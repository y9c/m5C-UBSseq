#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright © 2021 Ye Chang yech1990@gmail.com
# Distributed under terms of the MIT license.
#
# Created: 2021-10-06 01:53

"""pre-filter sites.

- V3: select sites with both unique and multi mapped reads.
      change pandas into polars to speed up.
- V4: read combined file directly.
"""

import argparse
import sys

import polars as pl

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

TOTAL_DEPTH = 20
TOTAL_SUPPORT = 3
AVERAGE_UNC_RATIO = 0.02
AVERAGE_CLU_RATIO = 0.5
AVERAGE_MUL_RATIO = 0.2


dfs = []
for f in args.input_files:
    df = (
        pl.read_ipc(f)
        .filter(
            (pl.col("d") >= TOTAL_DEPTH)
            & (pl.col("u") >= TOTAL_SUPPORT)
            & (pl.col("ur") >= AVERAGE_UNC_RATIO)
            & (pl.col("cr") < AVERAGE_CLU_RATIO)
            & (pl.col("mr") < AVERAGE_MUL_RATIO)
        )
        .select(["ref", "pos", "strand"])
    )
    print(f"Read data for {f}...", file=sys.stderr)
    dfs.append(df)

pl.concat(dfs, how="vertical").unique(maintain_order=True).write_csv(
    args.output_file, separator="\t", include_header=False
)
