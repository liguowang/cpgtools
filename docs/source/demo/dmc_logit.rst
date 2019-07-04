dmc_logit.py
=============

Description
-----------
This program performs differential CpG analysis using `logistic regression model <https://en.wikipedia.org/wiki/Logistic_regression>`_
based on proportion values. It allows for covariable analysis. Users can choose to use
"binomial" or "quasibinomial" family to model the data. The quasibinomial family estimates 
an addition parameter indicating the amount of the oversidpersion.

Options
------------

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
                        Sample IDs shoud match to the "Data file".
  -f FAMILY_FUNC, --family-FAMILY_FUNC
                        Error distribution and link function to be used in the
                        GLM model. Can be integer 1 or 2 with 1 -
                        "quasibinomial" and 2 - "binomial". Default-1.
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file.

Input files (examples)
------------------------

- `test_04_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_04_TwoGroup.tsv.gz>`_
- `test_04_TwoGroup.grp.csv <https://sourceforge.net/projects/cpgtools/files/test/test_04_TwoGroup.grp.csv>`_

Command
----------
::

 $ dmc_logit.py -i test_04_TwoGroup.tsv.gz -g test_04_TwoGroup.grp.csv -o output_quasibin
 $ dmc_logit.py -i test_04_TwoGroup.tsv.gz -g test_04_TwoGroup.grp.csv -f 2  -o output_bin

