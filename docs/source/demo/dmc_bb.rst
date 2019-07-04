dmc_bb.py
===========

Description
------------

This program performs differential CpG analysis using "beta binomial" model on proportion
values. It allows for covariant analysis. 

**Notes**
- You must install R package `aod <https://cran.r-project.org/web/packages/aod/index.html>`_ before running this program.


Options
--------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file-INPUT_FILE
                        Data file containing methylation proportions
                        (represented by "methyl_count,total_count", eg.
                        "20,30") with the 1st row containing sample IDs (must
                        be unique) and the 1st column containing CpG positions
                        or probe IDs (must be unique). This file can be a
                        regular text file or compressed file (.gz, .bz2).
  -g GROUP_FILE, --group-GROUP_FILE
                        Group file defining the biological groups of each
                        sample as well as other covariables such as gender,
                        age. The first varialbe is grouping variable (must be
                        categorical), all the other variables are considered
                        as covariates (can be categorial or continuous).
                        Sample IDs shoud match to the "Data file"..
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file.
                         
Input files
--------------

- `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
- `test_05_TwoGroup.grp.csv <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp.csv>`_

Command
------------
::

 $ python3 ../bin/dmc_bb.py -i test_04_TwoGroup.tsv.gz -g test_04_TwoGroup.grp.csv -o OUT_bb                       
