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
  -i INPUT_FILE, --input_file=INPUT_FILE
                        Data file containing methylation proportions
                        (represented by "methyl_count,total_count", eg.
                        "20,30") with the 1st row containing sample IDs (must
                        be unique) and the 1st column containing CpG positions
                        or probe IDs (must be unique). This file can be a
                        regular text file or compressed file (*.gz, *.bz2) or
                        accessible url.
  -g GROUP_FILE, --group=GROUP_FILE
                        Group file defining the biological groups of each
                        sample. It is a comma-separated 2 columns file with
                        the 1st column containing sample IDs, and the 2nd
                        column containing group IDs.  It must have a header
                        row. Sample IDs should match to the "Data file".
  -o OUT_FILE, --output=OUT_FILE
                        The prefix of the output file.
                        
Input files (examples)
-----------------------

- `test_09.tsv.gz <https://sourceforge.net/projects/cpgtools/files/test/test_09.tsv.gz/download>`_
- `test_09.grp.csv <https://sourceforge.net/projects/cpgtools/files/test/test_09.grp.csv/download>`_

Commands
---------
::

 $ dmc_fisher.py -i test_09.tsv.gz -g test_09.grp.csv -o test_fisher


Output
---------

- 3 columns ("Odds ratio", "pvalue" and "FDR adjusted pvalue") will append to the original
  table.
::

 $ head -5 test_fisher.pval.txt
 ID	LTS_MCR-1008	LTS_MCR-1035	STS_MCR-1021	STS_MCR-1251	OddsRatio	pval	adj.pval
 chr10:100011340	12,14	26,37	0,18	10,24	9.353846153846154	1.2116597355208375e-06	6.343768248800197e-05
 chr10:100011341	0,21	0,54	0,26	0,19	nan	1.0	1.0
 chr10:100011387	0,14	0,40	0,20	0,24	nan	1.0	1.0
 chr10:100011388	18,18	47,54	19,23	18,19	1.2548262548262548	0.7574366471769988	1.0
 chr10:100026933	16,30	28,55	7,40	13,19	2.0926829268292684	0.04119183894184185	0.2617016451197068
  
