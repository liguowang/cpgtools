dmc_nonparametric.py
======================

Description
-----------
This program performs differential CpG analysis uisng the  `Mann-Whitney U test <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html>`_
for two group comparison, and the `Kruskal-Wallis H-test <https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance>`_
for multiple groups comparison.

Options
-----------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file-INPUT_FILE
                        Data file containing beta values with the 1st row
                        containing sample IDs (must be unique) and the 1st
                        column containing CpG positions or probe IDs (must be
                        unique). Except for the 1st row and 1st column, any
                        non-numerical values will be considered as "missing
                        values" and ignored. This file can be a regular text
                        file or compressed file (.gz, .bz2).
  -g GROUP_FILE, --group-GROUP_FILE
                        Group file defining the biological group of each
                        sample. It is a comma-separated 2 columns file with
                        the 1st column containing sample IDs, and the 2nd
                        column containing group IDs. It must have a header
                        row. Sample IDs should match to the "Data file". Note:
                        automatically switch to use  Kruskal-Wallis H-test if
                        more than 2 groups were defined in this file.
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file.
                        

Input files (examples)
------------------------

- `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_
- `test_05_TwoGroup.grp.csv <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.grp.csv>`_
- `test_06_ThreeGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_06_ThreeGroup.tsv.gz>`_
- `test_06_ThreeGroup.grp.csv <https://sourceforge.net/projects/cpgtools/files/test/test_06_ThreeGroup.grp.csv>`_

Command
-----------
::
 
 $dmc_nonparametric.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv -o U_test
 
 $dmc_nonparametric.py -i test_06_TwoGroup.tsv.gz -g test_06_TwoGroup.grp.csv -o H_test

