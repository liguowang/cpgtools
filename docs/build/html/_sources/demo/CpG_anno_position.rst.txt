CpG_anno_position.py
=====================

Description
------------

This program adds annotation information to each CpG based on its **genomic position**.

Notes
------

- Input CpG (-i) and annotation (-a) BED files must have at least three columns, and must based on the same genome assembly   version
- If multiple regions from the annotation BED file are overlapped with the **same**
  CpG site, their names will be concatenated together.
- Since the input (-i) is a regular BED foramt file, this module can be uesd to annotate any genomic regions of interest. 

Pre-computed datasets
----------------------

`hg19_ENCODE_338TF_130Cell_E3.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_338TF_130Cell_E3.bed.gz>`_ (File size  = 108.2 MB)
	Transcription factor (TF) binding sites identified from `ChIP-seq <https://en.wikipedia.org/wiki/ChIP-sequencing>`_ experiments performed
	by the `ENCODE <https://www.encodeproject.org/>`_ project. Peaks from 1264 experiments
	representing 338 transcription factors in 130 cell types are combined (N = 10,560,472). BED format file was
	downloaded from the `UCSC Tabel Browser <http://genome.ucsc.edu/cgi-bin/hgTables>`_.
	
`hg19_ENCODE_DNaseI_125Cells_V3.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_DNaseI_125Cells_V3.bed.gz>`_ (File size  = 24.3 MB)
	`DNase I hypersensitivity sites <https://en.wikipedia.org/wiki/DNase_I_hypersensitive_site>`_ identified from
	`ENCODE <https://www.encodeproject.org/>`_ `DNase-seq <https://en.wikipedia.org/wiki/DNase-Seq>`_ experiments. Peaks
	from 125 cell types are combined (N = 1,867,665). BED format file was downloaded from
	the `UCSC Tabel Browser <http://genome.ucsc.edu/cgi-bin/hgTables>`_.
	
`hg19_ENCODE_chromHMM_states_9Cells.merge.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_chromHMM_states_9Cells.merge.bed.gz>`_ (File size  = 32.7 MB)
	Chromatin State Segmentation by `chromHMM <https://www.nature.com/articles/nmeth.1906>`_ 
	from `ENCODE <https://www.encodeproject.org/>`_. Chromatin states across 9 cell
	types (GM12878, H1-hESC, K562, HepG2, HUVEC, HMEC, HSMM, NHEK, NHLF) were learned by
	integrating 9 factors (CTCF, H3K27ac, H3K27me3, H3K36me3, H3K4me1,
	H3K4me2, H3K4me3, H3K9ac, H4K20me1 ) plus input. A total of 15 states were identified,
	include: State-1 (Active Promoter), state-2 (Weak Promoter), state-3 (Inactive/poised
	Promoter), state-4 and 5 (Strong enhancer), state-6 and 7 (Weak/poised enhancer),
	state-8 (insulator), state-9 (Transcriptional transition), state-10 (Transcriptional
	elongation), state-11 (Weak transcribed), state-12 (Polycomb-repressed), state-13
	(Heterochromatin or low signal), state-14 and 15 (Repetitive/Copy Number Variation).
	The Original chromatin state BED file was downloaded from the `UCSC Tabel Browser <http://genome.ucsc.edu/cgi-bin/hgTables>`_.

`hg19_FANTOM_enhancers_phase_1_and_2.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_FANTOM_enhancers_phase_1_and_2.bed.gz>`_
	`PHANTOM5 <http://fantom.gsc.riken.jp/5/>`_ human permissive enhancers downloaded from `here <http://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_and_2_expression_tpm_matrix.txt.gz>`_.

`hg19_ENCODE_H3K4me1_11_cellLines_ChIP.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_H3K4me1_11_cellLines_ChIP.bed.gz>`_ (File size  = 12.2 MB)
	H3K4me1 (marker of active and primed enhancer) peaks identified from `ENCODE <https://www.encodeproject.org/>`_ histone ChIP-seq experiments. Peaks from 11 cell
	types (GM12878, H1-hESC, HMEC, HSMM, HUVEC, HeLaS3, HepG2, K562, Monocytes-CD14+_RO01746,
	NHEK, NHLF) are combined (N = 1,435,550)

`hg19_ENCODE_H3K4me3_11_cellLines_ChIP.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_H3K4me3_11_cellLines_ChIP.bed.gz>`_ (File size  = 4.5 MB)
	H3K4me3 (marker of promoter) peaks identified from `ENCODE <https://www.encodeproject.org/>`_ histone ChIP-seq experiments. Peaks from 11 cell
	types (GM12878, H1-hESC, HMEC, HSMM, HUVEC, HeLaS3, HepG2, K562, Monocytes-CD14+_RO01746,
	NHEK, NHLF) are combined (N = 525,824)

`hg19_ENCODE_H3K27ac_11_cellLines_ChIP.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_H3K27ac_11_cellLines_ChIP.bed.gz>`_ (File size  = 5.7 MB)
	H3K27ac (marker of active enhancer) peaks identified from `ENCODE <https://www.encodeproject.org/>`_ histone ChIP-seq experiments. Peaks from 11 cell
	types (GM12878, H1-hESC, HMEC, HSMM, HUVEC, HeLaS3, HepG2, K562, Monocytes-CD14+_RO01746,
	NHEK, NHLF) are combined (N = 665,650)


These BED files were lifted over to hg38/GRCh38 using `CrossMap <http://crossmap.sourceforge.net/>`_.
The hg38-based annotation files are available from `here <https://sourceforge.net/projects/cpgtools/files/data/>`_ 

Options
-------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file=INPUT_FILE
                        Input CpG file in BED3+ format.
  -a ANNO_FILE, --annotation=ANNO_FILE
                        Input annotation file in BED3+ format.
  -w WINDOW_SIZE, --window=WINDOW_SIZE
                        Size of window centering on the middle-point of each
                        genomic region defined in the annotation BED file
                        (i.e., window_size*0.5 will be extended to up- and
                        down-stream from the middle point of each genomic
                        region). default=100
  -o OUT_FILE, --output=OUT_FILE
                        The prefix of the output file.
  -l, --header          If True, the first row of input CpG file is header.
                        default=False 


Input files (examples)
----------------------

- `test_01.hg19.bed6 <https://sourceforge.net/projects/cpgtools/files/test/test_01.hg19.bed6>`_
- `hg19_ENCODE_338TF_130Cell_E3.bed.gz <https://sourceforge.net/projects/cpgtools/files/data/hg19_ENCODE_338TF_130Cell_E3.bed.gz/download>`_

Command
-------

::
 
 
 $CpG_anno_position.py  -l -a hg19_ENCODE_338TF_130Cell_E3.bed.gz -i test_01.hg19.bed6 -o output
 

Output files
-------------

- output.anno.txt

::

 $ head -5 output.anno.txt
 #Chrom	Start	End	Name	Beta	Strand	hg19_ENCODE_338TF_130Cell_E3.bed
 chr1	10847	10848	cg26928153	0.8965	+	N/A
 chr1	10849	10850	cg16269199	0.7915	+	N/A
 chr1	15864	15865	cg13869341	0.9325	+	N/A
 chr1	534241	534242	cg24669183	0.7941	+	FOXA2,MNT
 
