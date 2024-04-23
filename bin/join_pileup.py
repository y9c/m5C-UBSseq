#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2024 Ye Chang yech1990@gmail.com
# Distributed under terms of the GNU license.
#
# Created: 2024-02-09 13:41


import polars as pl


def import_df(file_name, suffix):
    df = pl.read_csv(
        file_name,
        separator="\t",
        columns=["ref", "pos", "strand", "convertedBaseCount", "unconvertedBaseCount"],
        dtypes={
            "ref": pl.Utf8,
            "pos": pl.Int64,
            "strand": pl.Utf8,
            "convertedBaseCount": pl.Int64,
            "unconvertedBaseCount": pl.Int64,
        },
    )
    df = df.rename(
        {
            "convertedBaseCount": "convertedBaseCount_" + suffix,
            "unconvertedBaseCount": "unconvertedBaseCount_" + suffix,
        }
    )
    return df


def combine_files(*files):
    suffixes = [
        "unfiltered_uniq",
        "unfiltered_multi",
        "filtered_uniq",
        "filtered_multi",
    ]
    f = files[0]
    s = suffixes[0]
    df_com = import_df(f, s)
    for f, s in zip(files[1:], suffixes[1:]):
        df = import_df(f, s)
        df_com = df_com.join(df, on=["ref", "pos", "strand"], how="outer_coalesce")
    return df_com.fill_null(0)


if __name__ == "__main__":
    import argparse

    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-i",
        "--input-files",
        nargs=4,
        required=True,
        help="4 input files: unfiltered_uniq, unfiltered_multi, filtered_uniq, filtered_multi",
    )
    arg_parser.add_argument("-o", "--output-file", help="output file")
    args = arg_parser.parse_args()

    # Write the combined DataFrame to a CSV file
    combine_files(*args.input_files).write_ipc(args.output_file, compression="lz4")
