#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-02-09 13:41


import polars as pl


def import_df(file_name, suffix):
    count_cols = [
        "convertedBaseCount_unfiltered_uniq",
        "unconvertedBaseCount_unfiltered_uniq",
        "convertedBaseCount_unfiltered_multi",
        "unconvertedBaseCount_unfiltered_multi",
        "convertedBaseCount_filtered_uniq",
        "unconvertedBaseCount_filtered_uniq",
        "convertedBaseCount_filtered_multi",
        "unconvertedBaseCount_filtered_multi",
    ]
    df = pl.read_ipc(file_name).rename({col: col + "_" + suffix for col in count_cols})
    return df


def combine_files(*files):
    samples = [f.split("/")[-1].rsplit(".", 2)[0] for f in files]
    f = files[0]
    s = samples[0]
    df_com = import_df(f, s)
    for f, s in zip(files[1:], samples[1:]):
        df = import_df(f, s)
        df_com = df_com.join(df, on=["ref", "pos", "strand"], how="outer_coalesce")

    df_com = (
        df_com.with_columns(
            u=pl.sum_horizontal(
                f"unconvertedBaseCount_filtered_uniq_{s}" for s in samples
            ),
            d=pl.sum_horizontal(
                f"{t}_filtered_uniq_{s}"
                for s in samples
                for t in ["convertedBaseCount", "unconvertedBaseCount"]
            ),
            _t=pl.sum_horizontal(
                f"{t1}_unfiltered_{t2}_{s}"
                for s in samples
                for t1 in ["convertedBaseCount", "unconvertedBaseCount"]
                for t2 in ["uniq", "multi"]
            ),
        )
        .with_columns(
            # ur: unconverted ratio
            ur=pl.col("u") / pl.col("d"),
            # mr: multiple mapping ratio
            mr=pl.sum_horizontal(
                f"{t}_unfiltered_multi_{s}"
                for s in samples
                for t in ["convertedBaseCount", "unconvertedBaseCount"]
            )
            / pl.col("_t"),
            # cr: cluster ratio
            cr=1
            - pl.sum_horizontal(
                f"{t1}_filtered_{t2}_{s}"
                for s in samples
                for t1 in ["convertedBaseCount", "unconvertedBaseCount"]
                for t2 in ["uniq", "multi"]
            )
            / pl.col("_t"),
        )
        .with_columns(
            [
                pl.col(f"unconvertedBaseCount_filtered_uniq_{s}").alias(f"u_{s}")
                for s in samples
            ]
            + [
                (
                    pl.col(f"unconvertedBaseCount_filtered_uniq_{s}")
                    + pl.col(f"convertedBaseCount_filtered_uniq_{s}")
                ).alias(f"d_{s}")
                for s in samples
            ]
        )
        .select(
            ["ref", "pos", "strand", "u", "d", "ur", "mr", "cr"]
            + [f"{t}_{s}" for s in samples for t in ["u", "d"]]
        )
        .fill_null(0)
    )

    return df_com


if __name__ == "__main__":
    import argparse

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-i",
        "--input-files",
        nargs="+",
        required=True,
        help="Multiple input files to be combined",
    )
    arg_parser.add_argument("-o", "--output-file", help="output file")
    args = arg_parser.parse_args()

    # Write the combined DataFrame to a CSV file
    combine_files(*args.input_files).write_ipc(args.output_file, compression="lz4")
