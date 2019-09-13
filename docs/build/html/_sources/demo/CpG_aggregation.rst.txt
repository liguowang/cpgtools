CpG_aggregation.py
===================

Description
-------------

Aggregate proportion values of a list of CpGs that located in give genomic regions
(eg. CpG islands, promoters, exons, etc.).

**Example of input file**
::

 Chrom	Start	End	score
 chr1	100017748	100017749	3,10
 chr1	100017769	100017770	0,10
 chr1	100017853	100017854	16,21

**Notes**

Outlier CpG will be removed if the probability of observing its proportion value is less
than p-cutoff. For example, if alpha set to 0.05, and there are 10 CpGs (n - 10) located in a
particular genomic region, the p-cutoff of this genomic region is 0.005 (0.05/10). Supposing
the total reads mapped to this region is 100, out of which 25 are methylated reads (i.e.
regional methylation level (beta) - 25/100 - 0.25)

- The probability of observing CpG (3,10) is : `pbinom(q=3, size=10, prob=0.25) = 0.7759`
- The probability of observing CpG (0,10) is : `qpbinom(q=0, size=10, prob=0.25) = 0.05631`
- The probability of observing CpG (16,21) is : `pbinom(q=16, size=21, prob=0.25, lower.tail=F) = 1.19e-07` (outlier)

Options
-------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input=INPUT_FILE
                        Input CpG file in BED format. The first 3 columns
                        contain "Chrom", "Start", and "End". The 4th column
                        contains proportion values.
  -a ALPHA_CUT, --alpha=ALPHA_CUT
                        The chance of mistakingly assign a particular CpG as
                        an outlier for each genomic region. default-0.05
  -b BED_FILE, --bed=BED_FILE
                        BED3+ file specifying the genomic regions.
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.

Input files (examples)
----------------------

- `test_03_RRBS.bed.gz <https://sourceforge.net/projects/cpgtools/files/test/test_03_RRBS.bed.gz>`_
- `hg19.RefSeq.union.1Kpromoter.bed.gz <https://sourceforge.net/projects/cpgtools/files/test/hg19.RefSeq.union.1Kpromoter.bed.gz>`_

Command
--------

::

 $CpG_aggregation.py -b hg19.RefSeq.union.1Kpromoter.bed.gz  -i 0_du145_133_glp_sh1.bed -o out

Output
-------


::

 chr1    567292  568293  3       0       93      3       0       93
 chr1    713567  714568  6       0       100     6       0       100
 chr1    762401  763402  7       0       110     7       0       110
 chr1    762470  763471  10      0       158     10      0       158
 chr1    854571  855572  2       12      16      2       12      16
 chr1    860620  861621  16      91      232     16      91      232
 chr1    894178  895179  12      151     229     41      506     735 

**Description**

- Column1-3: Genome coordinates

- Column4-6: numbers of "CpG", "aggregated methyl reads", and "aggregate total reads" **after**  outlier filtering

- Column7-9: numbers of "CpG", "aggregated methyl reads", and "aggregate total reads" **before**  outlier filtering
