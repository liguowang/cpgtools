dmc_ttest.py
============

Description
------------
Differential CpG analysis using `T test <https://en.wikipedia.org/wiki/Student%27s_t-test>`_
for two groups comparison or `ANOVA <https://en.wikipedia.org/wiki/Analysis_of_variance>`_ 
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
                        column containing group IDs.  It must have a header
                        row. Sample IDs should match to the "Data file". Note:
                        automatically switch to use ANOVA if more than 2
                        groups were defined in this file.
  -p, --paired          If '-p/--paired' flag was specified, use paired t-test
                        which requires the equal number of samples in both
                        groups. Paired sampels are matched by the order. This
                        option will be ignored for multiple group analysis.
  -w, --welch           If '-w/--welch' flag was specified, using Welch's
                        t-test which does not assume the two samples have
                        equal variance.  If omitted, use standard two-sample
                        t-test (i.e. assuming the two samples have equal
                        variance). This option will be ignored for paired
                        t-test and multiple group analysis.
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
 
 #Two group comparison. Compare normal livers to HCV-related cirrhosis livers 
 $dmc_ttest.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv -o ttest_2G
 
 #Three group comparison. Compare normal livers, HCV-related cirrhosis livers, and liver cancers 
 $dmc_ttest.py -i test_06_ThreeGroup.tsv.gz -g test_06_ThreeGroup.grp.csv -o ttest_3G
 
Output files
---------------

- ttest_2G.pval.txt
- ttest_3G.pval.txt

