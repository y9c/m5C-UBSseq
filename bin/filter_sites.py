#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-10-06 01:53


import argparse

import polars as pl
from scipy.stats import binomtest

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument("-i", "--input-file", help="Input site file")
arg_parser.add_argument("-m", "--mask-file", help="mask file")
arg_parser.add_argument("-b", "--background-file", help="background file")
arg_parser.add_argument("-o", "--output-file", help="output file")


args = arg_parser.parse_args()

df_site = (
    pl.read_ipc(args.input_file)
    .with_columns(
        u=pl.col("unconvertedBaseCount_filtered_uniq"),
        d=pl.col("convertedBaseCount_filtered_uniq")
        + pl.col("unconvertedBaseCount_filtered_uniq"),
    )
    .with_columns(ur=pl.col("u") / pl.col("d"))
)

df_pre = pl.read_csv(
    args.mask_file,
    separator="\t",
    has_header=False,
    new_columns=["ref", "pos", "strand"],
    dtypes={"ref": pl.Utf8, "pos": pl.Int64, "strand": pl.Utf8},
)

bg_ratio = (
    df_site.join(df_pre, on=["ref", "pos", "strand"], how="anti")
    .get_column("ur")
    .drop_nans()
    .mean()
)
with open(args.background_file, "w") as f:
    f.write(f"{bg_ratio}\n")


def testp(successes, trials, p):
    if successes == 0 or trials == 0:
        return 1.0
    return binomtest(successes, trials, p, alternative="greater").pvalue


df_filter = (
    df_pre.join(df_site, on=["ref", "pos", "strand"], how="left")
    .with_columns(pl.col("u").fill_null(strategy="zero"))
    .with_columns(pl.col("d").fill_null(strategy="zero"))
    .select(["ref", "pos", "strand", "u", "d", "ur"])
    .with_columns(
        pval=pl.struct(["u", "d"]).map_elements(
            lambda x: testp(x["u"], x["d"], bg_ratio)
        )
    )
    .with_columns(
        passed=(pl.col("pval") < 0.001)
        & (pl.col("u") >= 2)
        & (pl.col("d") >= 10)
        & (pl.col("ur") > 0.02)
    )
)


df_filter.write_csv(args.output_file, separator="\t", include_header=True)
