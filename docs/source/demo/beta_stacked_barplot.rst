beta_stacked_barplot.py
=========================

Description
------------
This program creates stacked barplot for each sample. The stacked barplot showing
the proportions of CpGs whose beta values are falling into these 4 ranges:
1. [0.00,  0.25]        #first quantile
2. [0.25,  0.50]        #second quantile
3. [0.50,  0.75]        #third quantile
4. [0.75,  1.00]        #forth quantile

**Example of input file**

::

 CpG_ID  Sample_01       Sample_02       Sample_03       Sample_04
 cg_001  0.831035        0.878022        0.794427        0.880911
 cg_002  0.249544        0.209949        0.234294        0.236680


Options
-----------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file-INPUT_FILE
                        Data frame file containing beta values with the 1st
                        row containing sample IDs and the 1st column
                        containing CpG IDs.
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file.

Input files (examples)
---------------------------

- `cirrHCV_vs_normal.data.tsv <https://sourceforge.net/projects/cpgtools/files/test/cirrHCV_vs_normal.data.tsv>`_
                        
Command
--------
::

 $beta_stacked_barplot.py -i cirrHCV_vs_normal.data.tsv -o stacked_bar
 
Output files
---------------

- stacked_bar.r
- stacked_bar.pdf

 
.. image:: ../_static/stacked_bar.png
   :height: 600 px
   :width: 650 px
   :scale: 100 %  

