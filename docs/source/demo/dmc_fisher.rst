dmc_fisher.py
=============

Description
------------
This program performs differential CpG analysis using Fisher exact test on proportion value.
It applies to two sample comparison with no biological/technical replicates. If biological/
technical replicates are provided, methyl reads and total reads of all replicates will be
merged (i.e. ignores biological/technical variations)

Input file format
--------------------
::

 # number before "," indicates number of methyl reads, and number after "," indicates
 # number of total reads
 cgID        sample_1    sample_2
 CpG_1       129,170     166,178
 CpG_2       24,77       67,99

Options
----------


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
                        Group file defining the biological group of each
                        sample. It is a comma-separated two columns file with
                        the 1st column containing sample IDs, and the 2nd
                        column containing group IDs.  It must have a header
                        row. Sample IDs should match to the "Data file".
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file.
                        

Output
---------

- 3 columns ("Odds ratio", "pvalue" and "FDR adjusted pvalue") will append to the original
  table.
