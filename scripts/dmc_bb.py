#!/usr/bin/env python3
#
# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: wangliguo78@gmail.com
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.

"""
Description
-----------
This program performs differential CpG analysis using beta-binomial model
based on methylation proportions (in the form of "c,n", where "c" indicates
"Number of reads with methylated C", and "n" indicates "Number of total
reads". Both c and n are non-negative integers and c <= n).

Example of input
----------------
Below example showing input data on 2 CpGs of 3 groups (A, B, and C)
with each group has 3 replicates:

    cgID   A_1      A_2      A_3   B_1   B_2    B_3    C_1    C_2    C_3
    CpG_1  129,170  166,178  7,9   1     6,16   10,10  10,15  11,15  16,22  20,36
    CpG_2  0,77     0,99     0,85  0,77  1,37   3,37   0,42   0,153  0,6

Notes
-----
1. It can handle covariants.
2. Input is proportion values, not beta values.
3. You must install R package "aod" before running this program.
   https://cran.r-project.org/web/packages/aod/index.html
"""

import os
import re
import subprocess
import sys

from optparse import OptionParser

from cpgmodule import ireader
from cpgmodule._version import __version__
from cpgmodule.utils import printlog, read_grp_file2


def main(argv=None):
    """Command-line entry point."""
    usage = "%prog [options]\n"
    parser = OptionParser(usage, version="%prog " + __version__)

    parser.add_option(
        "-i",
        "--input_file",
        action="store",
        type="string",
        dest="input_file",
        help=(
            'Data file containing methylation proportions (represented by '
            '"methyl_count,total_count", eg. "20,30") with the 1st row '
            "containing sample IDs (must be unique) and the 1st column "
            "containing CpG positions or probe IDs (must be unique). This "
            "file can be a regular text file or compressed file (.gz, .bz2)"
        ),
    )
    parser.add_option(
        "-g",
        "--group",
        action="store",
        type="string",
        dest="group_file",
        help=(
            "Group file defining the biological groups of each sample as "
            "well as other covariables such as gender, age. The first "
            "variable is grouping variable (must be categorical), all the "
            "other variables are considered as covariates (can be "
            'categorical or continuous). Sample IDs should match to the '
            '"Data file".'
        ),
    )
    parser.add_option(
        "-o",
        "--output",
        action="store",
        type="string",
        dest="out_file",
        help="The prefix of the output file.",
    )

    options, args = parser.parse_args(args=argv)

    print()

    if not options.input_file:
        print(__doc__)
        parser.print_help()
        return 101

    if not options.group_file:
        print(__doc__)
        parser.print_help()
        return 102

    if not options.out_file:
        print(__doc__)
        parser.print_help()
        return 103

    if not os.path.isfile(options.input_file):
        print(
            'Input data file "%s" does not exist\n' % options.input_file
        )
        return 104

    if not os.path.isfile(options.group_file):
        print(
            'Input group file "%s" does not exist\n' % options.input_file
        )
        return 105

    rout = open(options.out_file + ".r", "w")

    print('library("aod")', file=rout)

    printlog('Read group file "%s" ...' % options.group_file)

    samples, cv_names, cvs, v_types = read_grp_file2(options.group_file)

    for cv_name in cv_names:
        print("%s: %s" % (cv_name, v_types[cv_name]))
        for sample in samples:
            print("\t" + sample + "\t" + cvs[cv_name][sample])

    print(
        "bbr1 <- function (cgid, m,t,%s){"
        % ",".join(cv_names),
        file=rout,
    )
    print(
        "\tdat <- data.frame(m=m, t=t, %s)"
        % ",".join(["=".join(i) for i in zip(cv_names, cv_names)]),
        file=rout,
    )
    print(
        "\tfit <- betabin(cbind(m,t - m) ~ %s, ~1, "
        'link=c("logit"), data=na.omit(dat))'
        % "+".join(cv_names),
        file=rout,
    )
    print("\ttest <- summary(fit)", file=rout)
    print("\tcoefs <- test@Coef$Estimate", file=rout)
    print('\tpvals = test@Coef$"Pr(> |z|)"', file=rout)
    print(
        "\tif(max(pvals, na.rm=T)>1){pvals = pvals + NA}",
        file=rout,
    )
    print(
        "\tif(sum(m, na.rm=T) == 0){pvals = pvals + NA}",
        file=rout,
    )
    print("\tnames = row.names(test@Coef)", file=rout)
    print('\tnames = gsub("2","",names)', file=rout)
    print(
        '\twrite.table(file="%s",'
        "x=matrix(c(cgid, as.vector(coefs), as.vector(pvals)), nrow=1),"
        'quote=FALSE, row.names=FALSE, sep="\\t", '
        'col.names=c("ID",paste(names, "coef",sep="."), '
        'paste(names, "pval",sep=".")))'
        % (options.out_file + ".results.txt"),
        file=rout,
    )
    print("}", file=rout)
    print("\n", file=rout)

    print(
        "bbr2 <- function (cgid, m,t,%s){"
        % ",".join(cv_names),
        file=rout,
    )
    print(
        "\tdat <- data.frame(m=m, t=t, %s)"
        % ",".join(["=".join(i) for i in zip(cv_names, cv_names)]),
        file=rout,
    )
    print(
        "\tfit <- betabin(cbind(m,t - m) ~ %s, ~1, "
        'link=c("logit"), data=na.omit(dat))'
        % "+".join(cv_names),
        file=rout,
    )
    print("\ttest <- summary(fit)", file=rout)
    print("\tcoefs <- test@Coef$Estimate", file=rout)
    print('\tpvals = test@Coef$"Pr(> |z|)"', file=rout)
    print(
        "\tif(max(pvals, na.rm=T)>1){pvals = pvals + NA}",
        file=rout,
    )
    print(
        "\tif(sum(m, na.rm=T) == 0){pvals = pvals + NA}",
        file=rout,
    )
    print("\tnames = row.names(test@Coef)", file=rout)
    print('\tnames = gsub("2","",names)', file=rout)
    print(
        '\twrite.table(file="%s",'
        "x=matrix(c(cgid, as.vector(coefs), as.vector(pvals)), nrow=1), "
        'quote=FALSE, row.names=FALSE, sep="\\t", '
        "col.names=FALSE, append=TRUE)"
        % (options.out_file + ".results.txt"),
        file=rout,
    )
    print("}", file=rout)
    print("\n", file=rout)

    printlog('Processing file "%s" ...' % options.input_file)

    line_num = 0

    for line in ireader.reader(options.input_file):
        line_num += 1
        fields = line.split()

        if len(fields) == 0:
            continue

        if line_num == 1:
            sample_ids = fields[1:]

            for sample in samples:
                if sample not in sample_ids:
                    printlog(
                        'Cannot find sample ID "%s" from file "%s"'
                        % (sample, options.input_file)
                    )
                    rout.close()
                    return 3

            for cv_name in cv_names:
                if v_types[cv_name] == "continuous":
                    print(
                        cv_name
                        + " <- c(%s)"
                        % ",".join(
                            [str(cvs[cv_name][s]) for s in sample_ids]
                        ),
                        file=rout,
                    )
                elif v_types[cv_name] == "categorical":
                    print(
                        cv_name
                        + " <- as.factor(c(%s))"
                        % ",".join(
                            [str(cvs[cv_name][s]) for s in sample_ids]
                        ),
                        file=rout,
                    )
                else:
                    printlog("unknown vaiable type!")
                    rout.close()
                    return 1

            print("\n", file=rout)
            continue

        methyl_reads = []
        total_reads = []
        cg_id = fields[0]

        for value in fields[1:]:
            match = re.match(r"(\d+)\s*\,\s*(\d+)", value)

            if match is None:
                methyl_reads.append("NaN")
                total_reads.append("NaN")
                continue

            c = int(match.group(1))
            n = int(match.group(2))

            if n >= c and n > 0:
                methyl_reads.append(c)
                total_reads.append(n)
            else:
                printlog("Incorrect data format!")
                print(fields)
                rout.close()
                return 1

        if line_num == 2:
            print(
                'bbr1("%s", c(%s), c(%s), %s)'
                % (
                    cg_id,
                    ",".join([str(read) for read in methyl_reads]),
                    ",".join([str(read) for read in total_reads]),
                    ",".join(cv_names),
                ),
                file=rout,
            )
        else:
            print(
                'bbr2("%s", c(%s), c(%s), %s)'
                % (
                    cg_id,
                    ",".join([str(read) for read in methyl_reads]),
                    ",".join([str(read) for read in total_reads]),
                    ",".join(cv_names),
                ),
                file=rout,
            )

    rout.close()

    try:
        printlog(
            'Runing Rscript file "%s" ...' % (options.out_file + ".r")
        )
        subprocess.call(
            "Rscript %s 2>%s"
            % (
                options.out_file + ".r",
                options.out_file + ".warnings.txt",
            ),
            shell=True,
        )
    except Exception:
        print(
            'Error: cannot run Rscript: "%s"'
            % (options.out_file + ".r"),
            file=sys.stderr,
        )
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
