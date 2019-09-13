beta_stats.py
==============

Description
-------------
This program gives basic information on CpGs located in each genomic region. It adds 6
columns to the input BED file:

1. Number of CpGs detected in the genomic region
2. Min methylation level
3. Max methylation level
4. Average methylation level across all CpGs
5. Median methylation level across all CpGs
6. Standard deviation

Options
----------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file=INPUT_FILE
                        BED6+ file specifying the C position. This BED file
                        should have at least six columns (Chrom, ChromStart,
                        ChromeEnd, Name, Beta_value, Strand).  Note: the first
                        base in a chromosome is numbered 0. This file can be a
                        regular text file or compressed file (.gz, .bz2)
  -r REGION_FILE, --region=REGION_FILE
                        BED3+ file of genomic regions. This BED file should
                        have at least 3 columns (Chrom, ChromStart,
                        ChromeEnd).
  -o OUT_FILE, --output=OUT_FILE
                        The prefix of the output file.
 
                        

Input files (examples)
-------------------------

- `test_02.bed6.gz <https://sourceforge.net/projects/cpgtools/files/test/test_02.bed6.gz>`_
- `hg19.RefSeq.union.1Kpromoter.bed <https://sourceforge.net/projects/cpgtools/files/test/hg19.RefSeq.union.1Kpromoter.bed.gz/download>`_


Command
-----------
::

 $beta_stats.py -r hg19.RefSeq.union.1Kpromoter.bed.gz -i test_02.bed6.gz -o region_stats

Output files
---------------

- region_stats.txt

