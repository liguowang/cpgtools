dmc_glm.py
===========

Description
------------

This program performs differential CpG analysis using `generalized liner model 
<https://en.wikipedia.org/wiki/Generalized_linear_model>`_. It allows
for covariants analysis.

Options
--------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file-INPUT_FILE
                        Data file containing beta values with the 1st row
                        containing sample IDs (must be unique) and the 1st
                        column containing CpG positions or probe IDs (must be
                        unique). This file can be regular text file or
                        compressed file (.gz, .bz2).
  -g GROUP_FILE, --group-GROUP_FILE
                        Group file defining the biological groups of each
                        sample as well as other covariables such as gender,
                        age. The first varialbe is grouping variable (must be
                        categorical), all the other variables are considered
                        as covariates (can be categorial or continuous).
                        Sample IDs shoud match to the "Data file".
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file.

Input files (examples)
------------------------

- `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
- `test_05_TwoGroup.grp.csv <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp.csv>`_
- `test_05_TwoGroup.grp2.csv <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp2.csv>`_
                        
Command
--------
::

 $dmc_glm.py  -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv -o GLM_2G
 
 $dmc_glm.py  -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp2.csv -o GLM_2G
 
Outpu files
--------------

- GLM_2G.results.txt
- GLM_2G.r
- GLM_2G.pval.txt (final results)
