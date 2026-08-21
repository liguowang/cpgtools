#!/usr/bin/env python3
#
# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: "deep.omics.lab@gmail.com"
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.


import sys,os
import pandas as pd
import logging
import subprocess
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import rpy2.robjects.packages as rpackages
from rpy2.robjects import pandas2ri
from dmc.utils import plot_corr


def DunedinPACE_clock(beta_file, outfile, metafile=None, delimiter=None,
                      ff='pdf', na_percent=0.2, ovr=False):
    try:
        devtools = importr('devtools')
    except:
        # Use "utils" to install R's "devtools" if not being done.
        utils = rpackages.importr('utils')
        utils.install_packages('devtools')

    # use "devtools" to install "DunedinPACE" package
    # https://github.com/danbelsky/DunedinPACE
    try:
        DunedinPACE = rpackages.importr('DunedinPACE')
    except:
        devtools.install_github("danbelsky/DunedinPACE", build_vignettes = False, quiet=True)

    if outfile is not None:
        out_prefix = outfile
    else:
        out_prefix = 'DunedinPACE_out'
    age_out = out_prefix + '.DNAm_age.tsv'
    r_out = out_prefix + '.plots.R'

    if ff.lower() in ['pdf', 'png']:
        scatter_out = out_prefix + '.scatter_plot.' + ff.lower()
    else:
        logging.error("Does not suppor format: %s!" % ff)
        sys.exit(0)

    outfiles = [ age_out, r_out, scatter_out]

    if ovr is True:
        logging.warning(
            "Over write existing files with prefix: %s" % out_prefix)
        for tmp in outfiles:
            try:
                os.remove(tmp)
            except FileNotFoundError:
                pass
    else:
        for tmp in outfiles:
            if os.path.exists(tmp):
                logging.error(
                    ("%s exists! Use different prefix or specify "
                    "\"--overwrite\" to replace existing files." % tmp))
                sys.exit(0)

    logging.info("Read input file: \"%s\"" % beta_file)
    input_df = pd.read_csv(beta_file, sep=delimiter, index_col=0, engine="python")
    (n_cpg, n_sample) = input_df.shape
    logging.info(
        "Input file: \"%s\", Number of CpGs: %d, Number of samples: %d" %
        (beta_file, n_cpg, n_sample))

    # calculate DunedinPACE
    with (ro.default_converter + pandas2ri.converter).context():
        r_from_pd_df = ro.conversion.get_conversion().py2rpy(input_df)
    a = DunedinPACE.PACEProjector(betas=r_from_pd_df, proportionOfProbesRequired=1-na_percent)
    names = list(input_df.columns)
    values = list(a[0][:])
    output = pd.DataFrame(data={'DunedinPACE' : values}, index=names)

    if metafile is not None:
        logging.info("Read meta information file: \"%s\"" % metafile)
        meta_df = pd.read_csv(metafile, sep=None, index_col=0, engine='python')
        meta_df.index = meta_df.index.astype(str)

        logging.info("Combining meta information with predicted age")
        output = pd.concat([output, meta_df], axis=1)
    
        # generate scatter plot between c_age and d_age
        c_age = []
        for col_id in output.columns:
            if col_id.lower() == 'age':
                c_age = output[col_id]
                break
        d_age = output['DunedinPACE']

        if len(c_age) >= 2 and len(c_age) == len(d_age):
            logging.info(
                "Writing R script of scatter plot. Save to: %s" % r_out)
            plot_corr(c_age, d_age, outfile=scatter_out, rfile=r_out)

    # save predicted age
    logging.info("Save predicted DNAm age to: %s" % age_out)
    output.to_csv(age_out, sep="\t", index_label="Sample_ID")

    logging.info("Running R script: %s" % r_out)
    try:
        subprocess.call("Rscript " + r_out, shell=True)
    except subprocess.CalledProcessError as e:
        print("Cannot generate pdf file from " + r_out, file=sys.stderr)
        print(e.output, file=sys.stderr)
        pass
    return output

if __name__=='__main__':
    DunedinPACE_clock(beta_file='Test1_blood_N20_EPICv1_beta.tsv', outfile=None, metafile='Test1_blood_N20_EPICv1_info.tsv')
