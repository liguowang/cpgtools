beta_m_conversion.py
=====================


Description
------------
Convert Beta-value into M-value or vice vers

**Example of input (beta)**

 CpG_ID	Sample_01	Sample_02	Sample_03	Sample_04
 cg_001	0.831035	0.878022	0.794427	0.880911
 cg_002	0.249544	0.209949	0.234294	0.236680
 cg_003	0.845065	0.843957	0.840184	0.824286

Options
----------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file=INPUT_FILE
                        Tab-separated data frame file containing beta values
                        with the 1st row containing sample IDs and the 1st
                        column containing CpG IDs. This file can be a regular
                        text file or compressed file (.gz, .bz2).
  -d DATA_TYPE, --dtype=DATA_TYPE
                        Input data type either "Beta" or "M".
  -o OUT_FILE, --output=OUT_FILE
                        The output file.

Input file (example)
--------------------

- `test_08.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_08.tsv.gz/download>`_

Command
-------
::

 $ beta_m_conversion.py  -i test_08.tsv -d Beta -o test_08_M.tsv
 
Output
-------
::

 $ head -5 test_08_M.tsv
 cg_ID	TCGA-BC-A10Q	TCGA-BC-A10R	TCGA-BC-A10S	TCGA-BC-A10T	TCGA-BC-A10U
 cg00000029	-0.9127840676229807	-0.6635535075463712	-0.9389653708375745	-1.1786876012968779	-0.6217264255944122
 cg00000165	-2.4833534763405667	-2.3330364850204406	-2.858145170950326	-2.914508967160336	-2.3645896606652745
 cg00000236	2.478873972561897	3.0777336083377693	2.6760378499862143	3.04301970048709	2.7096166677505145
 cg00000289	0.9943771370790748	0.13339998728363872	0.5981994318909333	1.2402989291699527	1.432741941887314
