from snakemake.utils import min_version
from collections import defaultdict

min_version("7.0")


configfile: "config.yaml"
workdir: "workspace"


CUSTOMIZED_GENES = []
BATCH = "inhibitor-AzaC"
SPECIES = "GRCh38"
LIBRARY_STRATEGY = "INLINE"
STRANDNESS = "F"
REF = config[f"ref_{SPECIES}"]
REFTYPES = ["genes", "genome", "contamination"]
REFTYPES_CALL = ["genes", "genome"]




read_ids = ["R1", "R2"]
pairend_run_ids = []
sample2data = defaultdict(dict)
group2sample = defaultdict(list)
for s, v in config[f"samples"].items():
    # picked treated samples only
    # if "treated" in s:
    if v.get("treated", True):
        group2sample[v["group"]].append(s)
    for i, v2 in enumerate(v["data"], 1):
        r = f"run{i}"
        sample2data[s][r] = {k3: os.path.expanduser(v3) for k3, v3 in v2.items()}
        if len(v2) == 2:
            pairend_run_ids.append(s + "_" + r)


rule all:
    input:
        "report_qc/cutadapt_qc.html",
        "report_qc/cutadaptPE_qc.html" if len(pairend_run_ids) > 0 else [],
        "report_qc/report_falco_after.html",
        expand("stat_reads/trimming/{sample}.tsv", sample=sample2data.keys()),
        expand("stat_reads/mapping/{sample}.tsv", sample=sample2data.keys()),
        "count_reads/genome_single.count",
        expand(
            "calculated_rate/{ref}_{is_filtered}.tsv.gz",
            ref=REFTYPES_CALL,
            is_filtered=["filtered", "unfiltered"],
        ),
        expand("filtered_table/{ref}.tsv.gz", ref=REFTYPES_CALL),


# fastq-join do not work very well with adpter reads
# fastp is better
rule join_pairend_reads:
    input:
        lambda wildcards: sample2data[wildcards.sample][wildcards.rn].values(),
    output:
        join=temp("merged_reads/{sample}_{rn}.fq.gz"),
        un1=temp("merged_reads/{sample}_{rn}.un1.fq.gz"),
        un2=temp("merged_reads/{sample}_{rn}.un2.fq.gz"),
    params:
        html=temp("merged_reads/{sample}_{rn}.fastp.html"),
        json=temp("merged_reads/{sample}_{rn}.fastp.json"),
    threads: 12
    resources:
        mem_mb=48000,
    run:
        if len(input) == 2:
            shell(
                """
        fastp --thread {threads} --merge --correction --overlap_len_require 10 --overlap_diff_percent_limit 20 -i {input[0]} -I {input[1]} --merged_out {output.join} --out1 {output.un1} --out2 {output.un2} -h {params.html} -j {params.json}
        """
            )
        else:
            shell(
                """
        ln -sf {input[0]} {output[0]}
        touch {output[1]}
        touch {output[2]}
        """
            )


# cut SE


rule cutadapt_SE:
    input:
        "merged_reads/{sample}_{rn}.fq.gz",
    output:
        fastq_cut="cut_adapter_SE/{sample}_{rn}.fq.gz",
        fastq_tooshort="cut_adapter_SE/{sample}_{rn}.tooshort.fq.gz",
        fastq_untrimmed="cut_adapter_SE/{sample}_{rn}.untrimmed.fq.gz"
        if LIBRARY_STRATEGY in ["INLINE"]
        else [],
        report1="cut_adapter_SE/{sample}_{rn}.cutadapt.step1.report",
        report2="cut_adapter_SE/{sample}_{rn}.cutadapt.step2.report",
    params:
        # with inline barcode
        barcode3=config["adapter"]["barcode3"],
        adapter3=config["adapter"]["adapter3"],
    threads: 12
    resources:
        mem_mb=48000,
    run:
        if LIBRARY_STRATEGY == "INLINE":
            # INLINE, double ligation method
            # WARNING: deal with NEB SMALLRNA kit in this chunk, but I should not use this method
            # p5 - 5ntUMI - insert - 5ntUMI - 6nt inline barcode - p7
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.barcode3}{params.adapter3};e=0.15;o=6;anywhere;" \
            --untrimmed-output={output.fastq_untrimmed} \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 5 -u -5 \
            --rename='{{id}}_{{cut_prefix}}{{cut_suffix}} {{comment}}' \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort} \
            -o {output.fastq_cut} - >{output.report2}
        """
            )
        elif LIBRARY_STRATEGY == "TAKARAV3":
            # p5 - 14nttUMI - reverse insert - p7
            shell(
                """
        seqtk seq {input} | \
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.adapter3};anywhere;o=4;e=0.15" \
            -o - - 2>{output.report1} | \
        cutadapt -j {threads} \
            -u -14 \
            --rename='{{id}}_{{cut_suffix}} {{comment}}' \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort} \
            -o {output.fastq_cut} - >{output.report2}
        """
            )
        elif LIBRARY_STRATEGY == "SWIFT":
            # p5 - [might be 6bp of polyC] - reverse insert (cDNA) - adaptase tail (CCCCCC) - p7
            # 6nt of polyG in 5' of R1 might from random RT primer
            # adaptase tail can be as long as 15bp at the 5' of R2 of polyG)
            # no UMI, but try to use random polyC tail as UMI
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.adapter3};e=0.15;o=4;anywhere" \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 6 -u -15 \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort} \
            -o {output.fastq_cut} - >{output.report2}
        """
            )
        elif LIBRARY_STRATEGY == "STRANDED":
            # VAHTS Stranded mRNA-seq Library Prep Kit
            # KAPA Stranded mRNA-Seq Kit (KAPA)
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.adapter3};e=0.15;o=4;anywhere" \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 2 -u -6 \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort} \
            -o {output.fastq_cut} - >{output.report2}
        """
            )
        else:
            # Small RNA, double ligation method, without barcode
            # p5 - insert - p7
            # trim 2nt on both end to increase quality
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.adapter3};e=0.15;o=6;anywhere;" \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 2 -u -2 \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort} \
            -o {output.fastq_cut} - >{output.report2}
        """
            )


rule report_cutadapt_SE:
    input:
        lambda wildcards: [
            f"cut_adapter_SE/{s}_{r}.cutadapt.step2.report"
            for s, v in sample2data.items()
            for r in v.keys()
        ],
    output:
        "report_qc/cutadapt_qc.html",
    shell:
        "multiqc -f -m cutadapt -n {output} {input}"


# cut PE
# The libraries of TAKARA v3 KIT is too long, so there is too much unjoined reads.
# run in PE mode to save the reads


rule cutadapt_PE:
    input:
        "merged_reads/{sample}_{rn}.un1.fq.gz",
        "merged_reads/{sample}_{rn}.un2.fq.gz",
    output:
        fastq_cut1="cut_adapter_PE/{sample}_{rn}_R1.fq.gz",
        fastq_cut2="cut_adapter_PE/{sample}_{rn}_R2.fq.gz",
        fastq_tooshort1="cut_adapter_PE/{sample}_{rn}_R1.tooshort.fq.gz",
        fastq_tooshort2="cut_adapter_PE/{sample}_{rn}_R2.tooshort.fq.gz",
        fastq_untrimmed1="cut_adapter_PE/{sample}_{rn}_R1.untrimmed.fq.gz"
        if LIBRARY_STRATEGY == "INLINE"
        else [],
        fastq_untrimmed2="cut_adapter_PE/{sample}_{rn}_R2.untrimmed.fq.gz"
        if LIBRARY_STRATEGY == "INLINE"
        else [],
        report1="cut_adapter_PE/{sample}_{rn}.cutadapt.step1.report",
        report2="cut_adapter_PE/{sample}_{rn}.cutadapt.step2.report",
    params:
        # with inline barcode
        barcode3=config["adapter"]["barcode3"],
        adapter_r1=config["adapter"]["smallRNA_r1"]
        if LIBRARY_STRATEGY == "SMALLRNA"
        else config["adapter"]["truseq_r1"],
        adapter_r2=config["adapter"]["smallRNA_r2"]
        if LIBRARY_STRATEGY in ["INLINE", "SMALLRNA"]
        else config["adapter"]["truseq_r2"],
    threads: 12
    resources:
        mem_mb=48000,
    run:
        if LIBRARY_STRATEGY == "INLINE":
            # INLINE, double ligation method
            # WARNING: deal with NEB SMALLRNA kit in this chunk, but I should not use this method
            # p5 - 5ntUMI - insert - 5ntUMI - 6nt inline barcode - p7
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.barcode3}{params.adapter_r1};e=0.15;o=6;anywhere;" \
            -A "{params.adapter_r2};e=0.15;o=6;anywhere;" \
            --untrimmed-output={output.fastq_untrimmed1} --untrimmed-paired-output={output.fastq_untrimmed2} \
            --interleaved \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 5 -U 5 \
            --rename='{{id}}_{{r1.cut_prefix}}{{r2.cut_prefix}} {{comment}}' \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort1} --too-short-paired-output={output.fastq_tooshort2} \
            --interleaved \
            -o {output.fastq_cut1} -p {output.fastq_cut2} - >{output.report2}
        """
            )
        if LIBRARY_STRATEGY == "TAKARAV3":
            # p5 - 14ntUMI - reverse insert - p7
            # NOTE: do not need to cut interleaved (-U -14), because both reads do not overlap
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.adapter_r1};anywhere;o=4;e=0.15" \
            -A "{params.adapter_r2};anywhere;o=4;e=0.15" \
            --interleaved \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -U 14 \
            --rename='{{id}}_{{r2.cut_prefix}} {{comment}}' \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort1} --too-short-paired-output={output.fastq_tooshort2} \
            --interleaved \
            -o {output.fastq_cut1} -p {output.fastq_cut2} - >{output.report2}
        """
            )
        elif LIBRARY_STRATEGY == "SWIFT":
            # p5 - [might be 6bp of polyC] - reverse insert (cDNA) - adaptase tail (CCCCCC) - p7
            # 6nt of polyG in 5' of R1 might from random RT primer
            # adaptase tail can be as long as 15bp at the 5' of R2 of polyG)
            # no UMI, but try to use random polyC tail as UMI
            # NOTE: do not need to cut interleaved, because both reads do not overlap
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.adapter_r1};anywhere;o=4;e=0.15" \
            -A "{params.adapter_r2};anywhere;o=4;e=0.15" \
            --interleaved \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 6 -U 15 \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort1} --too-short-paired-output={output.fastq_tooshort2} \
            --interleaved \
            -o {output.fastq_cut1} -p {output.fastq_cut2} - >{output.report2}
        """
            )
        elif LIBRARY_STRATEGY == "STRANDED":
            # INLINE, double ligation method
            # WARNING: deal with NEB SMALLRNA kit in this chunk, but I should not use this method
            # p5 - 5ntUMI - insert - 5ntUMI - 6nt inline barcode - p7
            # NOTE: do not need to cut interleaved, because both reads do not overlap
            # cut 2nt for quality
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.barcode3}{params.adapter_r1};e=0.15;o=6;anywhere;" \
            -A "{params.adapter_r2};e=0.15;o=6;anywhere;" \
            --untrimmed-output={output.fastq_untrimmed1} --untrimmed-paired-output={output.fastq_untrimmed2} \
            --interleaved \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 2 -u -2 -U 5 -U -2 \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort1} --too-short-paired-output={output.fastq_tooshort2} \
            --interleaved \
            -o {output.fastq_cut1} -p {output.fastq_cut2} - >{output.report2}
        """
            )
        else:
            # Small RNA, double ligation method, without barcode
            # p5 - insert - p7
            # trim 2nt on both end to increase quality
            shell(
                """
        cutadapt -j {threads} \
            -n 2 \
            -a "{params.adapter_r1};e=0.15;o=6;anywhere;" \
            -A "{params.adapter_r2};e=0.15;o=6;anywhere;" \
            --interleaved \
            -o - {input} 2>{output.report1} | \
        cutadapt -j {threads} \
            -u 2 -u -2 -U 2 -U -2 \
            --max-n=0 \
            -q 15 \
            --nextseq-trim=15 \
            -m 20 \
            --too-short-output={output.fastq_tooshort1} --too-short-paired-output={output.fastq_tooshort2} \
            --interleaved \
            -o {output.fastq_cut1} -p {output.fastq_cut2} - >{output.report2}
        """
            )


rule report_cutadapt_PE:
    input:
        lambda wildcards: [
            f"cut_adapter_PE/{s}_{r}.cutadapt.step2.report"
            for s, v in sample2data.items()
            for r in v.keys()
        ],
    output:
        "report_qc/cutadaptPE_qc.html",
    shell:
        "multiqc -f -m cutadapt -n {output} {input}"


# trimmed part qc


rule falco_after:
    input:
        "cut_adapter_SE/{sample}_{rn}.fq.gz",
    output:
        html="quality_control/falco_after/{sample}_{rn}/fastqc_report.html",
        text="quality_control/falco_after/{sample}_{rn}/fastqc_data.txt",
        summary="quality_control/falco_after/{sample}_{rn}/summary.txt",
    params:
        "quality_control/falco_after/{sample}_{rn}",
    shell:
        "falco -o {params} {input}"


rule report_falco_after:
    input:
        lambda wildcards: [
            f"quality_control/falco_after/{s}_{r}/fastqc_data.txt"
            for s, v in sample2data.items()
            for r in v.keys()
        ],
    output:
        "report_qc/report_falco_after.html",
    threads: 2
    resources:
        mem_mb=8000,
    shell:
        "multiqc -f -m fastqc -n {output} {input}"


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
    resources:
        mem_mb=56000,
    shell:
        """
        cat {input} >{output.fa}
        rm -f `dirname {output.index}`/`basename {output.index} ".CT.1.ht2"`.*.ht2
        ~/tools/hisat2/hisat-3n-build -p 12 --base-change C,T {output.fa} {params.index}
        """


# hisat2-3N (SE mapping mode)


rule hisat2_3n_mapping_genes:
    input:
        "cut_adapter_SE/{sample}_{rn}.fq.gz",
        "prepared_genes/genes.3n.CT.1.ht2" if CUSTOMIZED_GENES else [],
    output:
        sam=temp("run_mapping_SE/{sample}_{rn}.genes.sam"),
        fq=temp("run_mapping_SE/{sample}_{rn}.genes.fq"),
        summary="run_mapping_SE/{sample}_{rn}.genes.summary",
    params:
        hisat3n=config["path"]["hisat3n"],
        index=REF["genes"]["hisat3n"]
        if not CUSTOMIZED_GENES
        else "prepared_genes/genes",
        mapping="--directional-mapping"
        if STRANDNESS == "F"
        else "--directional-mapping-reverse"
        if STRANDNESS == "R"
        else "",
    threads: 24
    resources:
        mem_mb=56000,
    shell:
        """
        export TMPDIR="/scratch/midway3/yec"
        {params.hisat3n} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input[0]} {params.mapping} --all --norc --base-change C,T --mp 8,2 --no-spliced-alignment --un {output.fq} -S {output.sam}
        """


rule hisat2_3n_mapping_genome:
    input:
        "run_mapping_SE/{sample}_{rn}.genes.fq",
    output:
        sam=temp("run_mapping_SE/{sample}_{rn}.genome.sam"),
        fq=temp("run_mapping_SE/{sample}_{rn}.genome.fq"),
        summary="run_mapping_SE/{sample}_{rn}.genome.summary",
    params:
        hisat3n=config["path"]["hisat3n"],
        index=REF["genome"]["hisat3n"],
        mapping="--directional-mapping"
        if STRANDNESS == "F"
        else "--directional-mapping-reverse"
        if STRANDNESS == "R"
        else "",
    threads: 24
    resources:
        mem_mb=56000,
    shell:
        """
        export TMPDIR="/scratch/midway3/yec"
        {params.hisat3n} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input[0]} {params.mapping} --base-change C,T --pen-noncansplice 20 --mp 4,1 --un {output.fq} -S {output.sam}
        """


rule hisat2_3n_mapping_contamination:
    input:
        "run_mapping_SE/{sample}_{rn}.genome.fq",
    output:
        sam=temp("run_mapping_SE/{sample}_{rn}.contamination.sam"),
        fqz="run_unmapped/{sample}_{rn}.contamination.fq.gz",
        summary="run_mapping_SE/{sample}_{rn}.contamination.summary",
    params:
        hisat3n=config["path"]["hisat3n"],
        index=config["ref_contamination"]["hisat3n"],
        mapping="--directional-mapping"
        if STRANDNESS == "F"
        else "--directional-mapping-reverse"
        if STRANDNESS == "R"
        else "",
    threads: 24
    resources:
        mem_mb=56000,
    shell:
        """
        export TMPDIR="/scratch/midway3/yec"
        {params.hisat3n} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -U {input[0]} {params.mapping} --base-change C,T --mp 8,2 --no-spliced-alignment --un-gz {output.fqz} -S {output.sam}
        """


# hisat2-3N (PE mapping mode)


rule hisat2_3n_mapping_genes_PE:
    input:
        "cut_adapter_PE/{sample}_{rn}_R1.fq.gz",
        "cut_adapter_PE/{sample}_{rn}_R2.fq.gz",
        "prepared_genes/genes.3n.CT.1.ht2" if CUSTOMIZED_GENES else [],
    output:
        sam=temp("run_mapping_PE/{sample}_{rn}.genes.sam"),
        fq1=temp("run_mapping_PE/{sample}_{rn}_R1.genes.fq"),
        fq2=temp("run_mapping_PE/{sample}_{rn}_R2.genes.fq"),
        summary="run_mapping_PE/{sample}_{rn}.genes.summary",
    params:
        un="run_mapping_PE/{sample}_{rn}_R%.genes.fq",
        hisat3n=config["path"]["hisat3n"],
        index=REF["genes"]["hisat3n"]
        if not CUSTOMIZED_GENES
        else "prepared_genes/genes",
        mapping="--directional-mapping"
        if STRANDNESS == "F"
        else "--directional-mapping-reverse"
        if STRANDNESS == "R"
        else "",
    threads: 24
    resources:
        mem_mb=56000,
    shell:
        """
        export TMPDIR="/scratch/midway3/yec"
        {params.hisat3n} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input[0]} -2 {input[1]} {params.mapping} --all --norc --base-change C,T --mp 8,2 --no-spliced-alignment --un-conc {params.un} -S {output.sam}
        """


rule hisat2_3n_mapping_genome_PE:
    input:
        "run_mapping_PE/{sample}_{rn}_R1.genes.fq",
        "run_mapping_PE/{sample}_{rn}_R2.genes.fq",
    output:
        sam=temp("run_mapping_PE/{sample}_{rn}.genome.sam"),
        fq1=temp("run_mapping_PE/{sample}_{rn}_R1.genome.fq"),
        fq2=temp("run_mapping_PE/{sample}_{rn}_R2.genome.fq"),
        summary="run_mapping_PE/{sample}_{rn}.genome.summary",
    params:
        un="run_mapping_PE/{sample}_{rn}_R%.genome.fq",
        hisat3n=config["path"]["hisat3n"],
        index=REF["genome"]["hisat3n"],
        mapping="--directional-mapping"
        if STRANDNESS == "F"
        else "--directional-mapping-reverse"
        if STRANDNESS == "R"
        else "",
    threads: 24
    resources:
        mem_mb=56000,
    shell:
        """
        export TMPDIR="/scratch/midway3/yec"
        {params.hisat3n} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input[0]} -2 {input[1]} {params.mapping} --base-change C,T --pen-noncansplice 20 --mp 4,1 --un-conc {params.un} -S {output.sam}
        """


rule hisat2_3n_mapping_contamination_PE:
    input:
        "run_mapping_PE/{sample}_{rn}_R1.genome.fq",
        "run_mapping_PE/{sample}_{rn}_R2.genome.fq",
    output:
        sam=temp("run_mapping_PE/{sample}_{rn}.contamination.sam"),
        fq1="run_unmapped_PE/{sample}_{rn}_R1.contamination.fq.gz",
        fq2="run_unmapped_PE/{sample}_{rn}_R2.contamination.fq.gz",
        summary="run_mapping_PE/{sample}_{rn}.contamination.summary",
    params:
        un="run_unmapped_PE/{sample}_{rn}_R%.contamination.fq.gz",
        hisat3n=config["path"]["hisat3n"],
        index=config["ref_contamination"]["hisat3n"],
        mapping="--directional-mapping"
        if STRANDNESS == "F"
        else "--directional-mapping-reverse"
        if STRANDNESS == "R"
        else "",
    threads: 24
    resources:
        mem_mb=56000,
    shell:
        """
        export TMPDIR="/scratch/midway3/yec"
        {params.hisat3n} --index {params.index} -p {threads} --summary-file {output.summary} --new-summary -q -1 {input[0]} -2 {input[1]} {params.mapping} --base-change C,T --mp 8,2 --no-spliced-alignment --un-conc-gz {params.un} -S {output.sam}
        """


# join SE and PE mapping and sort


rule hisat2_3n_sort:
    input:
        "run_mapping_SE/{sample}_{rn}.{ref}.sam",
    output:
        "run_mapping_SE/{sample}_{rn}.{ref}.bam",
    params:
        samtools=config["path"]["samtools"],
    threads: 8
    resources:
        mem_mb=40000,
    shell:
        """
        {params.samtools} view -@ {threads} -F4 -b {input} | {params.samtools} sort -@ {threads} --write-index -m 4G -O BAM -o {output} -
        """


rule hisat2_3n_sort_PE:
    input:
        "run_mapping_PE/{sample}_{rn}.{ref}.sam",
    output:
        "run_mapping_PE/{sample}_{rn}.{ref}.bam",
    params:
        samtools=config["path"]["samtools"],
    threads: 8
    resources:
        mem_mb=40000,
    shell:
        """
        {params.samtools} view -@ {threads} -F4 -b {input} | {params.samtools} sort -@ {threads} --write-index -m 4G -O BAM -o {output} -
        """


rule join_SE_and_PE_mapping:
    input:
        "run_mapping_SE/{sample}_{rn}.{ref}.bam",
        "run_mapping_PE/{sample}_{rn}.{ref}.bam",
    output:
        temp("run_mapping_join_SE_and_PE/{sample}_{rn}.{ref}.bam"),
    params:
        samtools=config["path"]["samtools"],
    threads: 8
    resources:
        mem_mb=40000,
    shell:
        "{params.samtools} merge -@ {threads} -o {output} {input}"


################################################################################

# combine mapping results (multi run)


rule combine_runs:
    input:
        lambda wildcards: [
            f"run_mapping_join_SE_and_PE/{wildcards.sample}_{r}.{wildcards.ref}.bam"
            if len(v2) > 1
            else f"run_mapping_SE/{wildcards.s}_{r}.{wildcards.ref}.bam"
            for r, v in sample2data[wildcards.sample].items()
        ],
    output:
        "combined_mapping/{sample}.{ref}.bam",
    params:
        path_samtools=config["path"]["samtools"],
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        "{params.path_samtools} merge -@ {threads} -o {output} {input}"


# --input-fmt-option 'filter=[NH]==1'


## stat reads


rule count_cutadapt_reads:
    input:
        lambda wildcards: [
            os.path.join(x, f"{s}_{r}.cutadapt.step{i}.report")
            for r, v in sample2data[wildcards.sample].items()
            for x in (
                ["cut_adapter_SE", "cut_adapter_PE"]
                if len(v) > 1
                else ["cut_adapter_SE"]
            )
            for i in [1, 2]
        ],
    output:
        "stat_reads/trimming/{sample}.tsv",
    params:
        py=os.path.join(config["src_dir"], "parse_cutadapt_report.py"),
    shell:
        """
        {params.py} {input} >{output}
        """


## stat mapping


rule stat_mapping_number:
    input:
        bam=lambda wildcards: [
            f"combined_mapping/{wildcards.sample}.{ref}.bam" for ref in REFTYPES
        ],
    output:
        tsv="stat_reads/mapping/{sample}.tsv",
    params:
        path_samtools=config["path"]["samtools"],
        refs=REFTYPES,
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        # {params.path_samtools} flagstats -@ {threads} -O tsv $file | awk -v ref="$ref" '{{FS="\\t";OFS="\\t"}}$3 == "mapped"{{t=$1}}$3 == "primary mapped"{{p=$1}}END{{print ref,p; if(t > p)print ref"_multi",t-p}}' >> {output}
        """
        paste <(echo {params.refs} |  tr " " "\\n") <(echo {input.bam} |  tr " " "\\n") | while read ref file; do
            {params.path_samtools} view -@ {threads} -F 3980 -c $file | awk -v ref="$ref" '{{FS="\\t";OFS="\\t"}}NR==1{{print ref,$1}}' >> {output}
        done
        """


## remove duplicates


rule hisat2_3n_dedup:
    input:
        bam="combined_mapping/{sample}.{ref}.bam",
    output:
        bam="dedup_mapping/{sample}.{ref}.bam",
        log="dedup_mapping/{sample}.{ref}.log",
    params:
        path_umicollapse=config["path"]["umicollapse"],
    threads: 2
    resources:
        mem_mb=40000,
    run:
        if LIBRARY_STRATEGY in ["INLINE", "TAKARAV3"]:
            shell(
                """
            module load java
            export TMPDIR="/scratch/midway3/yec"
            java -server -Xms8G -Xmx36G -Xss100M -Djava.io.tmpdir=/scratch/midway3/yec -jar {params.path_umicollapse} bam \
                --two-pass -i {input.bam} -o {output.bam}  >{output.log}
            """
            )
        else:
            shell(
                """
                ln -sfr {input.bam} {output.bam}
                touch {output.log}
            """
            )


rule dedup_index:
    input:
        bam="dedup_mapping/{sample}.{ref}.bam",
    output:
        bai="dedup_mapping/{sample}.{ref}.bam.bai",
    threads: 6
    resources:
        mem_mb=24000,
    shell:
        """
        samtools index -@ {threads} {input}
        """


## gene expression


rule count_genome_single:
    input:
        expand("dedup_mapping/{sample}.genome.bam", sample=sample2data.keys()),
    output:
        "count_reads/genome_single.count",
    params:
        gtf=REF["genome"]["gtf"],
        name="gene_name"
    threads: 24
    resources:
        mem_mb=48000,
    shell:
        "featureCounts -T {threads} -O --largestOverlap -t exon -g {params.name} -a {params.gtf} -o {output} {input}"


########

## call mutation


rule hisat2_3n_calling_unfiltered_unique:
    input:
        "dedup_mapping/{sample}.{ref}.bam",
    output:
        "called_unfiltered/{sample}.{ref}.tsv.gz",
    params:
        hisat3ntable=config["path"]["hisat3ntable"],
        samtools=config["path"]["samtools"],
        fa=lambda wildcards: REF[wildcards.ref]["fa"]
        if wildcards.ref != "genes" or not CUSTOMIZED_GENES
        else "prepared_genes/genes.fa",
    threads: 24
    resources:
        mem_mb=48000,
    shell:
        """
        samtools view -e "rlen<100000" -h {input} | {params.hisat3ntable} -p {threads} -u --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | bgzip -@ {threads} -c > {output}
        """


# Update: previous version only call unique mapping reads


rule hisat2_3n_calling_unfiltered_multi:
    input:
        "dedup_mapping/{sample}.{ref}.bam",
    output:
        "called_unfiltered_multi/{sample}.{ref}.tsv.gz",
    params:
        hisat3ntable=config["path"]["hisat3ntable"],
        samtools=config["path"]["samtools"],
        fa=lambda wildcards: REF[wildcards.ref]["fa"]
        if wildcards.ref != "genes" or not CUSTOMIZED_GENES
        else "prepared_genes/genes.fa",
    threads: 24
    resources:
        mem_mb=48000,
    shell:
        # Update: previous version only call unique mapping reads
        """
        samtools view -e "rlen<100000" -h {input} | {params.hisat3ntable} -p {threads} -m --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | bgzip -@ {threads} -c > {output}
        """


## filter cluster effect and call mutation


rule hisat2_3n_filtering:
    input:
        bam="dedup_mapping/{sample}.{ref}.bam",
    output:
        converted=temp("hisat_converted/{sample}.{ref}.bam"),
        # unconverted=temp("hisat_unconverted/{sample}.{ref}.bam"),
    params:
        samtools=config["path"]["samtools"],
    threads: 2
    resources:
        mem_mb=6000,
    shell:
        # {params.samtools} view -@ {threads} -e "[XM] * 20 > (qlen-sclen) || [Zf] > 3 || 3 * [Zf] > [Zf] + [Yf]" {input.bam} -O BAM -o {output.unconverted}
        """
        {params.samtools} view -@ {threads} -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" {input.bam} -O BAM -o {output.converted}
        """


rule hisat2_3n_calling_filtered_unqiue:
    input:
        "hisat_converted/{sample}.{ref}.bam",
    output:
        "called_filtered/{sample}.{ref}.tsv.gz",
    params:
        hisat3ntable=config["path"]["hisat3ntable"],
        samtools=config["path"]["samtools"],
        fa=lambda wildcards: REF[wildcards.ref]["fa"]
        if wildcards.ref != "genes" or not CUSTOMIZED_GENES
        else "prepared_genes/genes.fa",
    threads: 24
    resources:
        mem_mb=48000,
    shell:
        """
        samtools view -e "rlen<100000" -h {input} | {params.hisat3ntable} -p {threads} -u --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | bgzip -@ {threads} -c > {output}
        """


# Update: previous version only call unique mapping reads
rule hisat2_3n_calling_filtered_multi:
    input:
        "hisat_converted/{sample}.{ref}.bam",
    output:
        "called_filtered_multi/{sample}.{ref}.tsv.gz",
    params:
        hisat3ntable=config["path"]["hisat3ntable"],
        samtools=config["path"]["samtools"],
        fa=lambda wildcards: REF[wildcards.ref]["fa"]
        if wildcards.ref != "genes" or not CUSTOMIZED_GENES
        else "prepared_genes/genes.fa",
    threads: 24
    resources:
        mem_mb=48000,
    shell:
        """
        samtools view -e "rlen<100000" -h {input} | {params.hisat3ntable} -p {threads} -m --alignments - --ref {params.fa} --output-name /dev/stdout --base-change C,T | bgzip -@ {threads} -c > {output}
        """


## calculate methy rate


rule calculate_methylation_rate:
    input:
        expand(
            "called_{{is_filtered}}/{sample}.{{ref}}.tsv.gz",
            sample=sample2data.keys(),
        ),
    output:
        "calculated_rate/{ref}_{is_filtered}.tsv.gz",
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        """
        (
        echo -e "Sample\\tChrom\\tPos\\tStrand\\tUnconverted\\tDepth\\tRatio"
        for file in {input}; do
            sample=`basename $file | cut -d. -f1`
            zcat $file | awk -v samplename="$sample" 'BEGIN{{FS="\\t"; OFS="\\t"}} NR > 1 {{ print samplename,$1,$2,$3,$7,$7+$5,$7/($7+$5) }}'
        done
        ) | bgzip -@ {threads} -c > {output}
        """


## picked sites and merge table


rule select_called_sites:
    input:
        # in paries
        lambda wildcards: [
            f"{t}/{s}.{wildcards.ref}.tsv.gz"
            for s in group2sample[wildcards.group]
            for t in ["called_filtered", "called_filtered_multi"]
        ],
    output:
        "selected_sites/{group}_{ref}.tsv.gz",
    params:
        # AVERAGE_DEPTH = 10
        # AVERAGE_RATIO = 0.02
        # TOTAL_SUPPORT = 3
        py=os.path.join(config["src_dir"], "select_called_sites_v3.py"),
    threads: 4
    resources:
        mem_mb=160000,
    shell:
        """
        {params.py} {output} {input}
        """


rule combined_select_sites:
    input:
        expand(
            "selected_sites/{group}_{{ref}}.tsv.gz",
            group=group2sample.keys(),
        ),
    output:
        "combined_sites/{ref}.tsv.gz",
    params:
        py=os.path.join(config["src_dir"], "combined_selected_sites.py"),
    resources:
        mem_mb=4000,
    shell:
        """
        {params.py} {output} {input}
        """


# pick background sites


rule pick_background_sites:
    input:
        call="called_filtered/{sample}.{ref}.tsv.gz",
        site="combined_sites/{ref}.tsv.gz",
    output:
        "picked_background/{sample}.{ref}.tsv.gz",
    shell:
        """
        zcat {input.call} | \
        awk 'BEGIN{{FS="\\t";OFS="\\t";print "Chrom", "Pos", "Strand", "Unconverted", "Depth"}} $5+$7>5{{print $1,$2,$3,$7,$5+$7}}' | \
        grep -w -F -v -f <(zcat {input.site} | cut -f 1-3) |\
        bgzip -c > {output}
        """


rule sumup_background_sites_by_group:
    input:
        # in paries
        lambda wildcards: [
            f"picked_background/{s}.{wildcards.ref}.tsv.gz"
            for s in group2sample[wildcards.group]
        ],
    output:
        "sumup_background/{group}_{ref}.tsv.gz",
    params:
        py=os.path.join(config["src_dir"], "sumup_background_ratio.py"),
    threads: 4
    resources:
        mem_mb=16000,
    shell:
        """
        {params.py} {output} {input}
        """


rule combined_background_sites:
    input:
        expand(
            "sumup_background/{group}_{{ref}}.tsv.gz",
            group=group2sample.keys(),
        ),
    output:
        "combined_background_by_group/{ref}.tsv.gz",
    resources:
        mem_mb=4000,
    shell:
        """
        echo {input} |  tr " " "\\n" | while read file; do
            sample=`basename $file | cut -d. -f1`
            zcat $file | awk -v sample="$sample" 'BEGIN{{FS="\\t";OFS="\\t"}}NR>1{{a+=$4;b+=$5}}END{{print sample,a,b,a/b}}' | bgzip -c >> {output}
        done
        """


rule stat_sample_background:
    input:
        expand("picked_background/{sample}.{{ref}}.tsv.gz", sample=sample2data.keys()),
    output:
        "combined_background/{ref}.tsv.gz",
    resources:
        mem_mb=4000,
    shell:
        """
        echo {input} |  tr " " "\\n" | while read file; do
            sample=`basename $file | cut -d. -f1`
            zcat $file | awk -v sample="$sample" 'BEGIN{{FS="\\t";OFS="\\t"}}NR>1{{a+=$4;b+=$5}}END{{print sample,a,b,a/b}}' | bgzip -c >> {output}
        done
        """


# pick putative sites


rule annotate_combined_sites:
    input:
        "combined_sites/{ref}.tsv.gz",
    output:
        "annotated_table/{ref}.tsv.gz",
    # params:
    # py=os.path.join(config["src_dir"], "annotate_merged_table.py"),
    resources:
        mem_mb=10000,
    shell:
        # variant version >= 0.0.36
        """
        zcat {input} | variant-effect -r {SPECIES} -s -H -c 1,2,3 | gzip -c > {output}
        """


rule merge_raw_data_to_sites:
    input:
        sites="annotated_table/{ref}.tsv.gz",
        raw=lambda wildcards: [
            f"{t}/{s}.{wildcards.ref}.tsv.gz"
            for s in sample2data.keys()
            for t in ["called_filtered", "called_filtered_multi"]
        ],
    output:
        "merged_table/{ref}.tsv.gz",
    params:
        py=os.path.join(config["src_dir"], "join_raw_table_v2.py"),
    resources:
        mem_mb=40000,
    shell:
        """
        {params.py} {output} {input.sites} {input.raw}
        """


rule filter_sites:
    input:
        table="merged_table/{ref}.tsv.gz",
        bg="combined_background/{ref}.tsv.gz",
    output:
        "filtered_table/{ref}.tsv.gz",
    params:
        py=os.path.join(config["src_dir"], "filter_table_by_pvalue.py"),
    shell:
        """
        {params.py} -b {input.bg} -i {input.table} -o {output}
        """
