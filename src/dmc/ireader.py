#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: deep.omics.lab@gmail.com
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See LICENSE.txt in the project root for the full license text.

import bz2
import gzip
import sys

from subprocess import PIPE, Popen
from urllib.request import urlopen


def nopen(fname, mode="rb"):
    """Open local, compressed, remote, stdin/stdout, or pipe input/output."""
    if not isinstance(fname, str):
        return fname

    if fname.startswith("|"):
        process = Popen(
            fname[1:],
            stdout=PIPE,
            stdin=PIPE,
            shell=True,
        )
        if mode.startswith("r"):
            return process.stdout
        return process

    if fname == "-":
        return {
            "r": sys.stdin,
            "w": sys.stdout,
        }[mode[0]]

    if fname.endswith((".gz", ".Z", ".z")):
        return gzip.open(fname, mode)

    if fname.endswith((".bz", ".bz2", ".bzip2")):
        return bz2.open(fname, mode)

    if fname.startswith(("http://", "https://", "ftp://")):
        return urlopen(fname)

    return open(fname, mode)


def reader(fname):
    """Yield decoded, stripped lines from an input source."""
    for line in nopen(fname):
        yield line.decode("utf-8").strip().replace("\r", "")
	
