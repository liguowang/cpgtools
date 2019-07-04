CpG_distrb_chrom.py
====================

Description
--------------

This program calculates the distribution of CpG over chromosomes

Options
--------

  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILES, --input-files-INPUT_FILES
                        Input CpG file(s) in BED3+ format. Multiple BED files
                        should be separated by "," (eg: "-i
                        file_1.bed,file_2.bed,file_3.bed"). BED file can be a
                        regular text file or compressed file (.gz, .bz2). The
                        barplot figures will NOT be generated if you provide
                        more than 12 samples (bed files). [required]
  -n FILE_NAMES, --names-FILE_NAMES
                        Shorter and meaningful names to label samples. Should
                        be separated by "," and match CpG BED files in number.
                        If not provided, basenames of CpG BED files will be
                        used to label samples. [optional]
  -s CHROM_SIZE, --chrom-size-CHROM_SIZE
                        Chromosome size file. Tab or space separated text file
                        with 2 columns: the first column is chromosome
                        name/ID, the second column is chromosome size. This
                        file will determine: (1) which chromosomes are
                        included in the final barplots, so do NOT include
                        'unplaced', 'alternative' contigs in this file. (2)
                        The order of chromosomes in the final barplots.
                        [required]
  -o OUT_FILE, --output-OUT_FILE
                        Prefix of the output file. [required]                        

Input files (examples)
-----------------------

- `450K_probe.hg19.bed3.gz <https://sourceforge.net/projects/cpgtools/files/test/450K_probe.hg19.bed3.gz>`_
- `850K_probe.hg19.bed3.gz <https://sourceforge.net/projects/cpgtools/files/test/850K_probe.hg19.bed3.gz>`_
- `hg19.chrom.sizes <https://sourceforge.net/projects/cpgtools/files/refgene/hg19.chrom.sizes>`_

Command
-----------
::
 
 $ chrom_distribution.py -i 450K_probe.hg19.bed3.gz,850K_probe.hg19.bed3.gz -n 450K,850K \
   -s hg19.chrom.sizes -o chromDist
 
**Output files**

- chromDist.txt
- chromDist.r
- chromDist.CpG_total.pdf
- chromDist.CpG_percent.pdf
- chromDist.CpG_perMb.pdf

Total CpG count per chromsome 

.. image:: ../_static/chromDist.CpG_total.png
   :height: 200 px
   :width: 500 px
   :scale: 100 %  

CpG percent on each chromosome (normalized to total CpGs)    

.. image:: ../_static/chromDist.CpG_percent.png
   :height: 200 px
   :width: 500 px
   :scale: 100 %  

CpG per Mb (normalized to chromsome size)   
 
.. image:: ../_static/chromDist.CpG_perMb.png
   :height: 200 px
   :width: 500 px
   :scale: 100 %  
