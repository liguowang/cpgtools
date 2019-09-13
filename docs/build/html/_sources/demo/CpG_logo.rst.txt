CpG_logo.py
============

Description
------------

This program generates a DNA motif logo for a given set of CpGs. To answer the question of
"what is the genomic context for a given list of CpGs ?". This program first extracts
genomic sequences around C position, and then generate `motif matrices <https://en.wikipedia.org/wiki/Position_weight_matrix>`_
include:

- position frequency matrix (PFM)
- position probability matrix (PPM)
- position weight matrix (PWM)
- `MEME <http://meme-suite.org/doc/meme-format.html>`_ format matrix
- `Jaspar <http://jaspar.genereg.net/>`_ format matrix

It also generates motif logo using `weblogo <https://github.com/WebLogo/weblogo>`_

**Notes**

- input BED file must have strand information.

Options
--------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file=INPUT_FILE
                        BED file specifying the C position. This BED file
                        should have at least six columns (Chrom, ChromStart,
                        ChromeEnd, name, score, strand).  Note: Must provide
                        correct *strand* information. This file can be a
                        regular text file or compressed file (.gz, .bz2).
  -r GENOME_FILE, --refgenome=GENOME_FILE
                        Reference genome seqeunces in FASTA format. Must be
                        indexed using the samtools "faidx" command.
  -e EXTEND_SIZE, --extend=EXTEND_SIZE
                        Number of bases extended to up- and down-stream.
                        default=5 (bp)
  -n MOTIF_NAME, --name=MOTIF_NAME
                        Motif name. default=motif
  -o OUT_FILE, --output=OUT_FILE
                        The prefix of the output file.
                        
Input files (examples)
------------------------

- Human reference genome sequences in FASTA format: `hg19.fa.gz 
  <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz>`_ and `hg38.fa.gz
  <http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz>`_
- `450_CH.hg19.bed.gz <https://sourceforge.net/projects/cpgtools/files/test/450_CH.hg19.bed.gz>`_                       

Command
-----------
::

 $CpG_logo.py -i 450_CH.hg19.bed.gz -r hg19.fa -o 450_CH

Output files
--------------

- 450_CH.logo.fa
- 450_CH.logo.jaspar
- 450_CH.logo.meme
- 450_CH.logo.pfm
- 450_CH.logo.ppm
- 450_CH.logo.pwm
- 450_CH.logo.logo.pdf

.. image:: ../_static/450_CH.logo.png
   :height: 400 px
   :width: 600 px
   :scale: 100 %  
