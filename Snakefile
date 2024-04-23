from snakemake.utils import min_version
from collections import defaultdict
from pathlib import Path

min_version("8.0")


configfile: Path(workflow.basedir) / "config.yaml"


workdir: "workspace"


BIN = config["path"]
REF = config["reference"]

CUSTOMIZED_GENES = [os.path.expanduser(i) for i in config.get("customized_genes", [])]
WITH_UMI = config.get("library", "") in ["INLINE", "TAKARAV3"]
MARKDUP = config.get("markdup", False)


SAMPLE2DATA = defaultdict(dict)
SAMPLE2LIB = defaultdict(dict)
GROUP2SAMPLE = defaultdict(list)
for s, v in config["samples"].items():
    if v.get("treated", True):
        # set default group as sample name, if not specified
        GROUP2SAMPLE[v.get("group", s)].append(s)
    SAMPLE2LIB[s] = v.get("library", config.get("library", ""))
    for i, v2 in enumerate(v.get("data", []), 1):
        r = f"run{i}"
        SAMPLE2DATA[s][r] = {k3: os.path.expanduser(v3) for k3, v3 in v2.items()}


INTERNALDIR = Path("internal_files")
TEMPDIR = Path(".tmp")

if os.environ.get("TMPDIR") is None:
    os.environ["TMPDIR"] = str(TEMPDIR)


envvars:
    "TMPDIR",


rule all:
    input:
        expand("report_reads/mapped/{sample}.tsv", sample=SAMPLE2DATA.keys()),
        [
            (
                INTERNALDIR / f"discarded_reads/{sample}_{rn}_R1.unmapped.fq.gz"
                if len(v) == 1
                else [
                    INTERNALDIR / f"discarded_reads/{sample}_{rn}_R1.unmapped.fq.gz",
                    INTERNALDIR / f"discarded_reads/{sample}_{rn}_R2.unmapped.fq.gz",
                ]
            )
            for sample, v in SAMPLE2DATA.items()
            for rn, v2 in v.items()
        ],
        expand(
            "detected_sites/filtered/{sample}.{ref}.tsv",
            sample=SAMPLE2DATA.keys(),
            ref=["genes", "genome"],
        ),
        # expand("prefilter_sites_combined/{ref}.tsv", ref=["genes", "genome"]),


# cut adapter


rule cutadapt_SE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
    output:
        fastq_cut=temp(TEMPDIR / "cut_adapter_SE/{sample}_{rn}_R1.fq.gz"),
        fastq_tooshort=INTERNALDIR / "discarded_reads/{sample}_{rn}_R1.tooshort.fq.gz",
        fastq_untrimmed=INTERNALDIR / "discarded_reads/{sample}_{rn}_R1.untrimmed.fq.gz",
        report="report_reads/trimming/{sample}_{rn}.json",
    params:
        library=lambda wildcards: SAMPLE2LIB[wildcards.sample],
    threads: 20
    shell:
        """
        cutseq -t {threads} -A {params.library} -m 20 --trim-polyA --ensure-inline-barcode --auto-rc -o {output.fastq_cut} -s {output.fastq_tooshort} -u {output.fastq_untrimmed} --json-file {output.report} {input} 
        """


rule cutadapt_PE:
    input:
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R1", "/"),
        lambda wildcards: SAMPLE2DATA[wildcards.sample][wildcards.rn].get("R2", "/"),
    output:
        fastq_cut=[
            temp(TEMPDIR / "cut_adapter_PE/{sample}_{rn}_R1.fq.gz"),
            temp(TEMPDIR / "cut_adapter_PE/{sample}_{rn}_R2.fq.gz"),
        ],
        fastq_tooshort=[
            INTERNALDIR / "discarded_reads/{sample}_{rn}_R1.tooshort.fq.gz",
            INTERNALDIR / "discarded_reads/{sample}_{rn}_R2.tooshort.fq.gz",
        ],
        fastq_untrimmed=[
            INTERNALDIR / "discarded_reads/{sample}_{rn}_R1.untrimmed.fq.gz",
            INTERNALDIR / "discarded_reads/{sample}_{rn}_R2.untrimmed.fq.gz",
        ],
        report="report_reads/trimming/{sample}_{rn}.report",
    params:
        library=lambda wildcards: SAMPLE2LIB[wildcards.sample],
    threads: 20
    shell:
        """
        cutseq -t {threads} -A {params.library} -m 20 --trim-polyA --ensure-inline-barcode --auto-rc -o {output.fastq_cut} -s {output.fastq_tooshort} -u {output.fastq_untrimmed} --json-file {output.report} {input} 
        """


# prepare genes index
# premap to rRNA, tRNA and other small RNA
# If study virus, then also premap to virus genome


rule prepare_genes_index:
    input:
        CUSTOMIZED_GENES,
    output:
        fa="prepared_genes/genes.fa",
        index="prepared_genes/genes.3n.CT.1.ht2",
    params:
        index="prepared_genes/genes",
    threads: 12
    shell:
        """
        cat {input} >{output.fa}
        rm -f `dirname {output.index}`/`basename {output.index} ".CT.1.ht2"`.*.ht2
        ~/tools/hisat2/hisat-3n-build -p 12 --base-change C,T {output.fa} {params.index}
        """


rule build_gene_index:
    input:
        fa="prepared_genes/genes.fa",
    output:
        fai="prepared_genes/genes.fa.fai",
    shell:
        """
        samtools faidx {output.fa}
        """


rule generate_saf_gene:
    input:
        fai="prepared_genes/genes.fa.fai",
    output:
        saf="prepared_genes/genes.saf",
    shell:
        """
        awk 'BEGIN{{OFS="\\t"}}{{print $1,$1,0,$2,"+"}}' {input} > {output}
        """


################################
# Mapping (SE mapping mode)
################################


rule hisat2_3n_mapping_contamination_SE:
    input:
        TEMPDIR / "cut_adapter_SE/{sample}_{rn}_R1.fq.gz",
    output:
        mapped=temp(TEMPDIR / "mapping_unsorted_SE/{sample}_{rn}.contamination.bam"),
        unmapped=temp(TEMPDIR / "mapping_discarded_SE/{sample}_{rn}.contamination.bam"),
        summary="report_reads/mapping/{sample}_{rn}.contamination.summary",
    params:
        index=REF["contamination"]["hisat3n"],
    threads: 24
    shell:
        """
        {BIN[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input[0]} --directional-mapping --base-change C,T --mp 8,2 --no-spliced-alignment | \
            {BIN[samtools]} view -@ {threads} -e '!flag.unmap' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule hisat2_3n_mapping_genes_SE:
    input:
        TEMPDIR / "unmapped_internal_SE/{sample}_{rn}_R1.contamination.fq.gz",
        "prepared_genes/genes.3n.CT.1.ht2" if CUSTOMIZED_GENES else [],
    output:
        mapped=temp(TEMPDIR / "mapping_unsorted_SE/{sample}_{rn}.genes.bam"),
        unmapped=temp(TEMPDIR / "mapping_discarded_SE/{sample}_{rn}.genes.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=(
            REF["genes"]["hisat3n"] if not CUSTOMIZED_GENES else "prepared_genes/genes"
        ),
    threads: 24
    shell:
        """
        {BIN[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input[0]} --directional-mapping --all --norc --base-change C,T --mp 8,2 --no-spliced-alignment | \
            {BIN[samtools]} view -@ {threads} -e '!flag.unmap' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule hisat2_3n_mapping_genome_SE:
    input:
        TEMPDIR / "unmapped_internal_SE/{sample}_{rn}_R1.genes.fq.gz",
    output:
        mapped=temp(TEMPDIR / "mapping_unsorted_SE/{sample}_{rn}.genome.bam"),
        unmapped=temp(TEMPDIR / "mapping_discarded_SE/{sample}_{rn}.genome.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genome.summary",
    params:
        index=REF["genome"]["hisat3n"],
    threads: 24
    shell:
        """
        {BIN[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input[0]} --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 | \
            {BIN[samtools]} view -@ {threads} -e '!flag.unmap' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule extract_unmap_bam_internal_SE:
    input:
        TEMPDIR / "mapping_discarded_SE/{sample}_{rn}.{reftype}.bam",
    output:
        temp(TEMPDIR / "unmapped_internal_SE/{sample}_{rn}_R1.{reftype}.fq.gz"),
    threads: 4
    shell:
        """
        {BIN[samtools]} fastq -@ {threads} -0 {output} {input}
        """


################################
# Mapping (PE mapping mode)
################################


rule hisat2_3n_mapping_contamination_PE:
    input:
        TEMPDIR / "cut_adapter_PE/{sample}_{rn}_R1.fq.gz",
        TEMPDIR / "cut_adapter_PE/{sample}_{rn}_R2.fq.gz",
    output:
        mapped=temp(TEMPDIR / "mapping_unsorted_PE/{sample}_{rn}.contamination.bam"),
        unmapped=temp(TEMPDIR / "mapping_discarded_PE/{sample}_{rn}.contamination.bam"),
        summary="report_reads/mapping/{sample}_{rn}.contamination.summary",
    params:
        index=REF["contamination"]["hisat3n"],
    threads: 24
    shell:
        """
        {BIN[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input[0]} -2 {input[1]} --directional-mapping --base-change C,T --mp 8,2 --no-spliced-alignment | \
            {BIN[samtools]} view -@ {threads} -e 'flag.proper_pair && !flag.unmap && !flag.munmap' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule hisat2_3n_mapping_genes_PE:
    input:
        TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R1.contamination.fq.gz",
        TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R2.contamination.fq.gz",
        "prepared_genes/genes.3n.CT.1.ht2" if CUSTOMIZED_GENES else [],
    output:
        mapped=temp(TEMPDIR / "mapping_unsorted_PE/{sample}_{rn}.genes.bam"),
        unmapped=temp(TEMPDIR / "mapping_discarded_PE/{sample}_{rn}.genes.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genes.summary",
    params:
        index=(
            REF["genes"]["hisat3n"] if not CUSTOMIZED_GENES else "prepared_genes/genes"
        ),
    threads: 24
    shell:
        """
        {BIN[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input[0]} -2 {input[1]} --directional-mapping --all --norc --base-change C,T --mp 8,2 --no-spliced-alignment | \
            {BIN[samtools]} view -@ {threads} -e 'flag.proper_pair && !flag.unmap && !flag.munmap' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule hisat2_3n_mapping_genome_PE:
    input:
        TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R1.genes.fq.gz",
        TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R2.genes.fq.gz",
    output:
        mapped=temp(TEMPDIR / "mapping_unsorted_PE/{sample}_{rn}.genome.bam"),
        unmapped=temp(TEMPDIR / "mapping_discarded_PE/{sample}_{rn}.genome.bam"),
        summary="report_reads/mapping/{sample}_{rn}.genome.summary",
    params:
        index=REF["genome"]["hisat3n"],
    threads: 24
    shell:
        """
        {BIN[hisat3n]} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input[0]} -2 {input[1]} --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 | \
            {BIN[samtools]} view -@ {threads} -e 'flag.proper_pair && !flag.unmap && !flag.munmap' -O BAM -U {output.unmapped} -o {output.mapped}
        """


rule extract_unmap_bam_internal_PE:
    input:
        TEMPDIR / "mapping_discarded_PE/{sample}_{rn}.{reftype}.bam",
    output:
        r1=temp(TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R1.{reftype}.fq.gz"),
        r2=temp(TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R2.{reftype}.fq.gz"),
    # wildcard_constraints: reftype="contamination|genes",
    threads: 4
    shell:
        """
        {BIN[samtools]} fastq -@ {threads} -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null -n {input}
        """


# keep unmap reads


ruleorder: extract_unmap_bam_final_PE > extract_unmap_bam_final_SE


rule extract_unmap_bam_final_SE:
    input:
        r1=TEMPDIR / "unmapped_internal_SE/{sample}_{rn}_R1.genome.fq.gz",
    output:
        r1=INTERNALDIR / "discarded_reads/{sample}_{rn}_R1.unmapped.fq.gz",
    threads: 4
    shell:
        """
        mv {input.r1} {output.r1}
        """


rule extract_unmap_bam_final_PE:
    input:
        r1=TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R1.genome.fq.gz",
        r2=TEMPDIR / "unmapped_internal_PE/{sample}_{rn}_R2.genome.fq.gz",
    output:
        r1=INTERNALDIR / "discarded_reads/{sample}_{rn}_R1.unmapped.fq.gz",
        r2=INTERNALDIR / "discarded_reads/{sample}_{rn}_R2.unmapped.fq.gz",
    threads: 4
    shell:
        """
        mv {input.r1} {output.r1}
        mv {input.r2} {output.r2}
        """


# sort and keep bam file for each run (for speedup reanalyse)


rule hisat2_3n_sort:
    input:
        lambda wildcards: TEMPDIR
        / (
            "mapping_unsorted_SE/{sample}_{rn}.{ref}.bam"
            if len(SAMPLE2DATA[wildcards.sample][wildcards.rn]) == 1
            else "mapping_unsorted_PE/{sample}_{rn}.{ref}.bam"
        ),
    output:
        INTERNALDIR / "run_sorted/{sample}_{rn}.{ref}.bam",
    threads: 16
    shell:
        """
        {BIN[samtools]} sort -@ {threads} --write-index -m 3G -O BAM -o {output} {input}
        """


################################################################################
# combine mapping results (multi run)
################################################################################


rule combine_runs:
    input:
        lambda wildcards: [
            INTERNALDIR / f"run_sorted/{wildcards.sample}_{r}.{wildcards.ref}.bam"
            for r in SAMPLE2DATA[wildcards.sample]
        ],
    output:
        temp(TEMPDIR / "combined_mapping/{sample}.{ref}.bam"),
    params:
        path_samtools=config["path"]["samtools"],
    threads: 8
    shell:
        "{params.path_samtools} merge -@ {threads} -o {output} {input}"


rule stat_mapping_number:
    input:
        bam=lambda wildcards: [
            TEMPDIR / f"combined_mapping/{wildcards.sample}.{ref}.bam"
            for ref in ["contamination", "genes", "genome"]
        ],
    output:
        tsv="report_reads/mapped/{sample}.tsv",
    params:
        refs=["contamination", "genes", "genome"],
    threads: 4
    shell:
        """
        paste <(echo {params.refs} |  tr " " "\\n") <(echo {input.bam} |  tr " " "\\n") | while read ref file; do
            {BIN[samtools]} view -@ {threads} -F 3980 -c $file | awk -v ref="$ref" '{{FS="\\t";OFS="\\t"}}NR==1{{print ref,$1}}' >> {output}
        done
        """


rule dedup_mapping:
    input:
        bam=TEMPDIR / "combined_mapping/{sample}.{ref}.bam",
    output:
        bam=INTERNALDIR / "aligned_bam/{sample}.{ref}.bam",
        txt="report_reads/dedup/{sample}.{ref}.log",
    params:
        tmp=os.environ["TMPDIR"],
    threads: 20
    run:
        if WITH_UMI:
            shell(
                """
            /software/java-15.0.2-el8-x86_64/bin/java -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir={params.tmp} -jar {BIN[umicollapse]} bam \
                -t 2 -T {threads} --data naive --merge avgqual --two-pass -i {input.bam} -o {output.bam} >{output.txt}
            """
            )
        elif MARKDUP:
            shell(
                """
                ~/tools/jdk8u322-b06-jre/bin/java -Xmx36G -jar ~/tools/gatk-4.2.5.0/gatk-package-4.2.5.0-local.jar MarkDuplicates \
                    -I {input} -O {output.bam} -M {output.txt} \
                    --DUPLICATE_SCORING_STRATEGY SUM_OF_BASE_QUALITIES --REMOVE_DUPLICATES true --VALIDATION_STRINGENCY SILENT --TMP_DIR {params.tmp}
            """
            )
        else:
            shell(
                """
                cp {input.bam} {output.bam}
                touch {output.txt}
            """
            )


rule dedup_index:
    input:
        bam=INTERNALDIR / "aligned_bam/{sample}.{ref}.bam",
    output:
        bai=INTERNALDIR / "aligned_bam/{sample}.{ref}.bam.bai",
    threads: 6
    shell:
        """
        {BIN[samtools]} index -@ {threads} {input}
        """


####################################
# call mutation
####################################


rule hisat2_3n_calling_unfiltered_unique:
    input:
        INTERNALDIR / "aligned_bam/{sample}.{ref}.bam",
    output:
        temp(TEMPDIR / "unfiltered_unique/{sample}.{ref}.tsv.gz"),
    params:
        fa=lambda wildcards: (
            REF[wildcards.ref]["fa"]
            if wildcards.ref != "genes" or not CUSTOMIZED_GENES
            else "prepared_genes/genes.fa"
        ),
    threads: 16
    shell:
        # ref     pos     strand  convertedBaseQualities  convertedBaseCount      unconvertedBaseQualities        unconvertedBaseCount
        """
        {BIN[samtools]} view -e "rlen<100000" -h {input} | {BIN[hisat3ntable]} -p {threads} -u --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | bgzip -@ {threads} -c > {output}
        """


rule hisat2_3n_calling_unfiltered_multi:
    input:
        INTERNALDIR / "aligned_bam/{sample}.{ref}.bam",
    output:
        temp(TEMPDIR / "unfiltered_multi/{sample}.{ref}.tsv.gz"),
    params:
        fa=lambda wildcards: (
            REF[wildcards.ref]["fa"]
            if wildcards.ref != "genes" or not CUSTOMIZED_GENES
            else "prepared_genes/genes.fa"
        ),
    threads: 16
    shell:
        """
        {BIN[samtools]} view -e "rlen<100000" -h {input} | {BIN[hisat3ntable]} -p {threads} -m --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | bgzip -@ {threads} -c > {output}
        """


## filter cluster effect and call mutation


rule hisat2_3n_filtering:
    input:
        INTERNALDIR / "aligned_bam/{sample}.{ref}.bam",
    output:
        temp(TEMPDIR / "hisat_converted/{sample}.{ref}.bam"),
    threads: 4
    shell:
        """
        {BIN[samtools]} view -@ {threads} -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" {input} -O BAM -o {output}
        """


rule hisat2_3n_calling_filtered_unqiue:
    input:
        TEMPDIR / "hisat_converted/{sample}.{ref}.bam",
    output:
        temp(TEMPDIR / "filtered_unique/{sample}.{ref}.tsv.gz"),
    params:
        fa=lambda wildcards: (
            REF[wildcards.ref]["fa"]
            if wildcards.ref != "genes" or not CUSTOMIZED_GENES
            else "prepared_genes/genes.fa"
        ),
    threads: 16
    shell:
        """
        {BIN[samtools]} view -e "rlen<100000" -h {input} | {BIN[hisat3ntable]} -p {threads} -u --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | bgzip -@ {threads} -c > {output}
        """


rule hisat2_3n_calling_filtered_multi:
    input:
        TEMPDIR / "hisat_converted/{sample}.{ref}.bam",
    output:
        temp(TEMPDIR / "filtered_multi/{sample}.{ref}.tsv.gz"),
    params:
        fa=lambda wildcards: (
            REF[wildcards.ref]["fa"]
            if wildcards.ref != "genes" or not CUSTOMIZED_GENES
            else "prepared_genes/genes.fa"
        ),
    threads: 16
    shell:
        """
        {BIN[samtools]} view -e "rlen<100000" -h {input} | {BIN[hisat3ntable]} -p {threads} -m --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | cut -f 1,2,3,5,7 | bgzip -@ {threads} -c > {output}
        """


rule join_pileup:
    input:
        lambda wildcards: [
            TEMPDIR / f"{t}/{wildcards.sample}.{wildcards.ref}.tsv.gz"
            for t in [
                "unfiltered_unique",
                "unfiltered_multi",
                "filtered_unique",
                "filtered_multi",
            ]
        ],
    output:
        INTERNALDIR / "count_sites/{sample}.{ref}.arrow",
    threads: 6
    shell:
        # "ref", "pos", "strand",
        # "convertedBaseCount_unfiltered_uniq",
        # "unconvertedBaseCount_unfiltered_uniq",
        # "convertedBaseCount_unfiltered_multi",
        # "unconvertedBaseCount_unfiltered_multi",
        # "convertedBaseCount_filtered_uniq",
        # "unconvertedBaseCount_filtered_uniq",
        # "convertedBaseCount_filtered_multi",
        # "unconvertedBaseCount_filtered_multi",
        """
        {BIN[join_pileup.py]} -i {input} -o {output}
        """


rule group_pileup:
    input:
        lambda wildcards: [
            INTERNALDIR / f"count_sites/{sample}.{wildcards.ref}.arrow"
            for sample in GROUP2SAMPLE[wildcards.group]
        ],
    output:
        INTERNALDIR / "group_sites/{group}.{ref}.arrow",
    threads: 6
    shell:
        # u: total unconverted count of unique filtered
        # d: total depth of unique filtered
        # ur: average unconverted ratio of unique filtered
        # mr: average multiple mapped ratio of unfiltered
        # cr: average cluster ratio of all
        # 
        # u_sample1: ...
        # d_sample1: ...
        # u_sample2: ...
        # d_sample2: ...
        # ...
        """
        {BIN[group_pileup.py]} -i {input} -o {output}
        """


# prefilter sites and merge table


rule combined_select_sites:
    input:
        expand(
            INTERNALDIR / "group_sites/{group}.{{ref}}.arrow",
            group=GROUP2SAMPLE.keys(),
        ),
    output:
        "detected_sites/prefilter/{ref}.tsv",
    shell:
        """
        {BIN[select_sites.py]} -i {input} -o {output}
        """


# pick putative sites by sample


rule stat_sample_background:
    input:
        site=INTERNALDIR / "count_sites/{sample}.{ref}.arrow",
        mask="detected_sites/prefilter/{ref}.tsv",
    output:
        background="detected_sites/background/{sample}.{ref}.tsv",
        filtered="detected_sites/filtered/{sample}.{ref}.tsv",
    threads: 2
    shell:
        """
        {BIN[filter_sites.py]} -i {input.site} -m {input.mask} -b {output.background} -o {output.filtered}
        """
