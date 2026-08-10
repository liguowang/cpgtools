#!/usr/bin/env python3

# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: wangliguo78@gmail.com
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.

"""Convert DNA methylation beta values to M values, or M values to beta values.

Features
--------
* Supports plain text and compressed input (.gz, .bz2, .xz).
* Preserves the CpG identifier column.
* Uses an epsilon offset to avoid infinite M values when beta values are 0 or 1.
* Writes tab-delimited output.

Conversion formulas
-------------------
Beta -> M

    M = log2(beta / (1 - beta))

M -> Beta

    Beta = 2**M / (2**M + 1)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from cpgmodule._version import __version__
from cpgmodule.utils import printlog


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=None, engine="python", compression="infer")


def main():
    parser = argparse.ArgumentParser(
        description="Convert DNA methylation beta values and M values.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i","--input_file",required=True,type=Path,
                        help="Input beta or M-value matrix.")
    parser.add_argument("-d","--dtype",choices=("Beta","M","beta","m"),
                        default="Beta",
                        help="Input matrix type.")
    parser.add_argument("-o","--out_prefix",required=True,
                        help="Output file prefix.")
    parser.add_argument("-e","--epsilon",type=float,default=1e-6,
                        help="Clamp beta values away from 0 and 1.")
    args = parser.parse_args()

    df = read_table(args.input_file)
    if df.shape[1] < 2:
        parser.error("Input must contain one CpG column and at least one sample column.")

    out = Path(args.out_prefix)

    ids = df.iloc[:,0]
    vals = df.iloc[:,1:].apply(pd.to_numeric, errors="coerce")

    if args.dtype.lower() == "beta":
        printlog(f'Converting Beta values in "{args.input_file}" to M values...')
        vals = np.log2(np.clip(vals, args.epsilon, 1-args.epsilon) /
                       (1 - np.clip(vals, args.epsilon, 1-args.epsilon)))
        outfile = out.with_suffix(".m.tsv")
    else:
        printlog(f'Converting M values in "{args.input_file}" to Beta values...')
        vals = (2.0 ** vals) / ((2.0 ** vals) + 1.0)
        outfile = out.with_suffix(".beta.tsv")

    result = pd.concat([ids, vals], axis=1)
    result.to_csv(outfile, sep="\t", index=False, float_format="%.6f")

    printlog(f'Saved converted matrix to "{outfile}".')


if __name__ == "__main__":
    raise SystemExit(main())
