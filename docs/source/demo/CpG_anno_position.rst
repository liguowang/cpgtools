CpG_anno_probe.py
==================

Description
-----------------

This program adds annotation information to each CpG based on its **genomic position**.

Notes
------

- Input CpG and BED files must have at least three columns, and must based on the same genome assembly version
- If multiple regions from the annotation BED file are overlapped with the **same**
  CpG site, their names will be concatenated together.
  
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
 
 # probe IDs are located in the 4th column (-p 3)
 
 $CpG_anno_probe.py -p 3 -l -a MethylationEPIC_CpGtools.tsv -i test_01.hg19.bed6 -o output
 
 or (take gzipped files as input) 
 
 $ CpG_anno_position.py -l -i test_01.hg19.bed6 -a hg19_ENCODE_338TF_130Cell_E3.bed.gz -o output



Output files
-------------

- output.anno.txt
