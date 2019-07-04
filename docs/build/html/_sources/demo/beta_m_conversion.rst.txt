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
  -i INPUT_FILE, --input-INPUT_FILE
                        Tab separated data frame file containing beta values
                        with the 1st row containing sample IDs and the 1st
                        column containing CpG IDs. This file can be a regular
                        text file or compressed file (.gz, .bz2) or
                        accessible url.
  -d DATA_TYPE, --dtype-DATA_TYPE
                        Input data type either "Beta" or "M".
  -o OUT_FILE, --output-OUT_FILE
                        Output file.
