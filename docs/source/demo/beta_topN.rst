beta_topN.py
=============

Description
------------
This program picks the top N rows (according to standard deviation) from the input file.
The resulting file can be used for clustering and PCA analysis.

**Example of input**

 CpG_ID  Sample_01       Sample_02       Sample_03       Sample_04
 cg_001  0.831035        0.878022        0.794427        0.880911
 cg_002  0.249544        0.209949        0.234294        0.236680
 cg_003  0.845065        0.843957        0.840184        0.824286

Options
-----------

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file=INPUT_FILE
                        Tab-separated data frame file containing beta values
                        with the 1st row containing sample IDs and the 1st
                        column containing CpG IDs.
  -c CPG_COUNT, --count=CPG_COUNT
                        Number of most variable CpGs (ranked by standard
                        deviation) to keep. default=1000
  -o OUT_FILE, --output=OUT_FILE
                        The prefix of the output file.

Input files (examples)
------------------------

- `test_05_TwoGroup.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_05_TwoGroup.tsv.gz>`_

Command
--------
::

 $beta_topN.py -i test_05_TwoGroup.tsv.gz -c 500 -o test_05_TwoGroup

Output file
------------

- test_05_TwoGroup.sortedStdev.tsv
- test_05_TwoGroup.sortedStdev.topN.tsv
