CpG_anno_probe.py
==================

Description
-----------------

This program adds comprehensive annotation information to each 450K/850K array probe ID.
It will add 17 columns to the original input data file. These 17 columns include
(from left to right):

+-----------------------+-------------------------------------------------------------------------+
| Header Name           |Description                                                              |
+-----------------------+-------------------------------------------------------------------------+
| hg19_pos              |The genomic position of the CpG on human genome assembly `hg19 (or       |
|                       |GRCh37) <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/>`_      |
+-----------------------+-------------------------------------------------------------------------+
| hg38_pos              |The genomic position of the CpG on human genome assembly `hg38 (or       |
|                       |GRCh38) <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/>`_.     |
+-----------------------+-------------------------------------------------------------------------+
| strand                |Strand of the CpG. Value - "R" (reverse strand) or "F" (forward strand). |
+-----------------------+-------------------------------------------------------------------------+
| geneSymbol            |Genes the CpG has been assigned to. "N/A" indicates no genes were found. |
|                       |This is retrieved from the Illumina `MethylationEpic v1.0 B4             |
|                       |<https://support.illumina.com/downloads/infinium-methylationepic-v1-0-   |
|                       |product-files.html>`_ manifest file.                                     |
+-----------------------+-------------------------------------------------------------------------+
| CpGisland             |The CpG island (CGI) that overlaps with this CpG. "N/A" indicates no     |
|                       |CGIs were found.                                                         |
+-----------------------+-------------------------------------------------------------------------+
| with_450K             |Boolean indicating whether this CpG probe is also included in 450K.      |
|                       |"0" - No, "1"- Yes.                                                      |
+-----------------------+-------------------------------------------------------------------------+
| SNP_ID                |SNPs (rsID) that are close to this CpG. Multiple SNPs are separated      |
|                       |by ";". "N/A" indicates no SNPs were found.                              |
+-----------------------+-------------------------------------------------------------------------+
| SNP_distance          |The nucleotide distances between SNPs and the CpG.                       |
+-----------------------+-------------------------------------------------------------------------+
| SNP_MAF               |The `minor allele frequencies (MAF) <https://en.wikipedia.org/wiki       |
|                       |/Minor_allele_frequency>`_ of SNPs.                                      |
+-----------------------+-------------------------------------------------------------------------+
| Cross_Reactive        |Boolean ("0" - No, "1"- Yes) indicating whether this CpG could be        |
|                       |affected by cross-hybridization or underlying genetic variation as       |
|                       |reported by this `paper <https://genomebiology.biomedcentral.com/        |
|                       |articles/10.1186/s13059-016-1066-1>`_.                                   |
+-----------------------+-------------------------------------------------------------------------+
| ENCODE_TF_ChIP        |Transcription factor (TF) binding sites identified from ChIP-seq         |
|                       |experiments performed by the `ENCODE <https://www.encodeproject.org/>`_  |
|                       |project. Peaks from 1264 experiments representing 338 transcription      |
|                       |factors in 130 cell types are combined (N - 10,560,472).                 |
|                       |BED format file was downloaded from the `UCSC Tabel Browser              |
|                       |<https://genome.ucsc.edu/cgi-bin/hgTables>`_, and a detailed description |
|                       |is provided `here <https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid-      |
|                       |732007223_QUJBO5BMeBu3R7xczOAWQ0UV9A1f&c-chr9&g-encRegTfbsClustered>`_.  |
+-----------------------+-------------------------------------------------------------------------+
| ENCODE_DNaseI         |DNase I hypersensitivity sites identified from ENCODE `DNase-seq         |
|                       |<https://en.wikipedia.org/wiki/DNase-Seq>`_ experiments. Peaks from      |
|                       |125 cell types are combined (N - 1,867,665). BED format file was         |
|                       |downloaded from the `UCSC Table Browser                                  |
|                       |<https://genome.ucsc.edu/cgi-bin/hgTables>`_, and a detailed description |
|                       |is provided `here <https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid-      |
|                       |732007223_QUJBO5BMeBu3R7xczOAWQ0UV9A1f&c-chr9&g-                         |
|                       |wgEncodeRegDnaseClustered>`_.                                            |
+-----------------------+-------------------------------------------------------------------------+
|ENCODE_H3K27ac_ChIP    |H3K27ac peaks identified from ENCODE histone ChIP-seq experiments. Peaks |
|                       |from 11 cell types (GM12878, H1-hESC, HMEC, HSMM, HUVEC, HeLaS3, HepG2,  |
|                       |K562, Monocytes-CD14+_RO01746, NHEK, NHLF) are combined (N - 665,650)    | 
+-----------------------+-------------------------------------------------------------------------+
|ENCODE_H3K4me1_ChIP    |H3K4me1 peaks identified from ENCODE histone ChIP-seq experiments. Peaks |
|                       |from 11 cell types (GM12878, H1-hESC, HMEC, HSMM, HUVEC, HeLaS3, HepG2,  |
|                       |K562, Monocytes-CD14+_RO01746, NHEK, NHLF) are combined (N - 1,435,550)  | 
+-----------------------+-------------------------------------------------------------------------+
|ENCODE_H3K4me3_ChIP    |H3K4me3 peaks identified from ENCODE histone ChIP-seq experiments. Peaks |
|                       |from 11 cell types (GM12878, H1-hESC, HMEC, HSMM, HUVEC, HeLaS3, HepG2,  |
|                       |K562, Monocytes-CD14+_RO01746, NHEK, NHLF) are combined (N - 525,824)    | 
+-----------------------+-------------------------------------------------------------------------+
|ENCODE_chromHMM        |Chromatin State Segmentation by `chromHMM <https://www.nature.com/       |
|                       |articles/nmeth.1906>`_ from ENCODE. Chromatin states across 9 cell types |
|                       |(GM12878,  H1-hESC, K562, HepG2, HUVEC, HMEC, HSMM, NHEK, NHLF) were     |
|                       |learned by computationally by integrating 9 factors (CTCF, H3K27ac,      |
|                       |H3K27me3, H3K36me3, H3K4me1, H3K4me2, H3K4me3, H3K9ac, H4K20me1 )        |
|                       |plus input. A total of 15 states were identified, include: State-1       |
|                       |(Active Promoter), state-2 (Weak Promoter), state-3 (Inactive/poised     |
|                       |Promoter), state-4 and 5 (Strong enhancer), state-6 and 7                |
|                       |(Weak/poised enhancer), state-8 (insulator), state-9 (Transcriptional    |
|                       |transition), state-10 (Transcriptional elongation), state-11 (Weak       |
|                       |transcribed), state-12 (Polycomb-repressed), state-13 (Heterochromatin or| 
|                       |low signal), state-14 and 15 (Repetitive/Copy Number Variation).         |
|                       |Orignal chromatin state BED file was downloaded from `UCSC Table Browser |
|                       |<https://genome.ucsc.edu/cgi-bin/hgTables>`_, and detailed description   |
|                       |is provided `here <https://genome.ucsc.edu/cgi-bin/hgTrackUi?hgsid-      |
|                       |732007223_QUJBO5BMeBu3R7xczOAWQ0UV9A1f&c-chr9&g-wgEncodeBroadHmm>`_.     |
+-----------------------+-------------------------------------------------------------------------+
|FANTOM_enhancer        |PHANTOM5 human enhancers downloaded from `here <http://fantom.gsc.riken. |
|                       |jp/5/datafiles/latest/extra/Enhancers/human_permissive_enhancers_phase_1_|
|                       |and_2_expression_tpm_matrix.txt.gz>`_.                                   |
+-----------------------+-------------------------------------------------------------------------+

Notes
-------

- For peaks identified from ENCODE ChIP-seq and DNase-seq (ENCODE_TF_ChIP, ENCODE_H3K27ac_ChIP,
  ENCODE_H3K4me1_ChIP, ENCODE_H3K4me3_ChIP, and ENCODE_DNaseI), we require the probe  must be
  located in the 100 bp window centered on the **middle** of the peak.

Options
-------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input_file=INPUT_FILE
                        Input data file (Tab-separated) with a certain column
                        containing 450K/850K array CpG IDs. This file can be
                        a regular text file or compressed file (.gz, .bz2).
  -a ANNO_FILE, --annotation=ANNO_FILE
                        Annotation file. This file can be a regular text file 
                        or compressed file (.gz, .bz2). 
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
  -p PROBE_COL, --probe_column=PROBE_COL
                        The number specifying which column contains probe IDs.
                        Note: the column index starts with 0. default-0.
  -l, --header          Input data file has a header row.
 


Input files (examples)
----------------------

- `test_01.hg19.bed6 <https://sourceforge.net/projects/cpgtools/files/test/test_01.hg19.bed6>`_
- `MethylationEPIC_CpGtools.tsv.gz <https://sourceforge.net/projects/cpgtools/files/data/MethylationEPIC_CpGtools.tsv.gz>`_

Command
-------

::
 
 # probe IDs are located in the 4th column (-p 3)
 
 $CpG_anno_probe.py -p 3 -l -a MethylationEPIC_CpGtools.tsv -i test_01.hg19.bed6 -o output
 
 or (take gzipped files as input) 
 
 $CpG_anno_probe.py -p 3 -l -a MethylationEPIC_CpGtools.tsv.gz -i test_01.hg19.bed6.gz -o output

 @ 2019-06-28 09:12:41: Read annotation file "../epic/MethylationEPIC_CpGtools.tsv" ...
 @ 2019-06-28 09:12:52: Add annotation information to "test_01.hg19.bed6" ... 

Output files
-------------

- output.anno.txt
