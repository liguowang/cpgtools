beta_selectNBest.py
===================

Description
------------
Select the K best features according to the K highest scores. Scores can be measured by:

* ANOVA F-value between label/feature for classification tasks.
* Mutual information for a discrete target.
* Chi-squared stats of non-negative features for classification tasks.

This is similar to `beta_topN.py <https://cpgtools.readthedocs.io/en/latest/demo/beta_topN.html>`_,
but requires a group file.

**Example of input**
::

 CpG_ID  Sample_01       Sample_02       Sample_03       Sample_04
 cg_001  0.831035        0.878022        0.794427        0.880911
 cg_002  0.249544        0.209949        0.234294        0.236680
 cg_003  0.845065        0.843957        0.840184        0.824286

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit
   -i INPUT_FILE, --input_file=INPUT_FILE
                         Tab-separated data frame file containing beta values
                         with the 1st row containing sample IDs and the 1st
                         column containing CpG IDs.
   -g GROUP_FILE, --group=GROUP_FILE
                         Comma-separated group file defining the biological
                         groups of each sample.
   -c CPG_COUNT, --topK=CPG_COUNT
                         Number of top features to select. default=100
   -s SCORE_FUNCTION, --score-function=SCORE_FUNCTION
                         Scoring function used to measure the dependency
                         between features scores and labels. Must be "chisq"
                         (chi-squared statistic), "anova" (ANOVA F-value), or
                         "mi" (mutual information). default=chisq
   -o OUT_FILE, --output=OUT_FILE
                         The prefix of the output file.

Input files (examples)
------------------------

- `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_

Command
--------
::

 $beta_selectNBest.py -i test_05_TwoGroup.tsv -g test_05_TwoGroup.grp.csv  -o test_05_out

 @ 2024-09-05 09:29:38: Reading input file: "test_05_TwoGroup.tsv"
 @ 2024-09-05 09:29:38: 0 rows with missing values were removed.
 @ 2024-09-05 09:29:38: Transposing data matrix ...
 @ 2024-09-05 09:29:38: Total number of features: 10000
 @ 2024-09-05 09:29:38: Reading group file: "test_05_TwoGroup.grp.csv"
 @ 2024-09-05 09:29:38: Using Chi Square statistic to select features ...
 @ 2024-09-05 09:29:39: Total number of selected features : 100
 @ 2024-09-05 09:29:39: Writing to file: "test_05_out.selectedFeatures.tsv"


Output file
------------

- test_05_out.selectedFeatures.tsv
