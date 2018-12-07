## CpGtools -- Tools to analyze and visualize DNA methylation data

## Part 1: Installation

### Part 1.1: Prerequisites
CpGtools are written in [Python](https://www.python.org/). In particular,  **Python3 (v3.5.x)**
is required to run all scripts in CpGtools. Some scripts also need **R** to generate graphs and 
run generalized linear model (GLM).  

- [Python 3](https://www.python.org/downloads/) and [pip3](https://pip.pypa.io/en/stable/installing/)
- [R](https://www.r-project.org/)

### Part 1.2: Python Dependencies
Note: these packages will be automatically installed when you use [pip3](https://pip.pypa.io/en/stable/installing/)
to install CpGtools.

- [numpy](http://www.numpy.org/)
- [scipy](https://www.scipy.org/)
- [pysam](https://pypi.org/project/pysam/)
- [bx-python](https://pypi.org/project/bx-python/)
- [pyBigWig](https://pypi.org/project/pyBigWig/)
- [sklearn](https://www.scilearn.com/)
- [weblogo](https://pypi.org/project/weblogo/)

### Part 1.3: Install [pip3](https://pip.pypa.io/en/stable/installing/) (Skip this step if you already have [pip3](https://pip.pypa.io/en/stable/installing/))

1. First, download **[get-pip.py](https://bootstrap.pypa.io/get-pip.py)**
		curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
2. Then run the following:
		python get-pip.py
3. Run the following code to check:
		pip3 --version
		pip 18.1 from /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/pip (python 3.6)
		
		which pip3
		/Library/Frameworks/Python.framework/Versions/3.6/bin/pip
		
		which pip	
		/Library/Frameworks/Python.framework/Versions/3.6/bin/pip
Note that **pip** is actualy a soft link to the same executable file path with **pip3**. so you can use pip directly. 

### Part 1.4: Install CpGtools

You can run the following command to install CpGtools and all its dependencies. 
	
	pip install cpgtools (not ready)		

You can run the following command to **upgrade** CpGtools and all its dependencies. 	
	
	pip install cpgtools --upgrade (not ready)

## Part 2: Usage Information

### annotate_CpG.py
---

#### Overview
This program annotate CpGs by assigning them to gene's regulatory domains. Follows the
"[Basel plus extension rules](http://great.stanford.edu/public/html/index.php)" used by GREAT:

**Basal regulatory domain**:
A gene's basal regulatory domain is a window around its TSS (transcription start site). In 
particular, basal regulatory domain is obtained by extending '-u' basepairs (default = 5 kb)
to the upstream and '-d' basepairs (default = 1 kb) to the downstream of TSS regardless of
other nearby genes.

**Extended regulatory domain**:
The gene's basal regulatory domain is further extended in both directions to the nearest gene's
basal regulatory domain but no more than the maximum extension (specified by '-e', default =
1000 kb) in one direction.	

#### Basic usage

```text

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        BED file specifying the C position. Must have at least
                        3 columns (chrom start end). Note: the first base in a
                        chromosome is numbered 0. BED file can be regular or
                        compressed by 'gzip' or 'bz'. [required]
  -r GENE_FILE, --refgene=GENE_FILE
                        Reference gene model in standard BED-12 format
                        (https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
  -u BASAL_UP_SIZE, --basal-up=BASAL_UP_SIZE
                        Size of extension to upstream of TSS (used to define
                        gene's "basal regulatory domain"). default=5000 (bp)
  -d BASAL_DOWN_SIZE, --basal-down=BASAL_DOWN_SIZE
                        Size of extension to downstream of TSS (used to define
                        gene's basal regulatory domain). default=1000 (bp)
  -e EXTENSION_SIZE, --extension=EXTENSION_SIZE
                        Size of extension to both up- and down-stream of TSS
                        (used to define gene's "extended regulatory domain").
                        default=1000000 (bp)
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of output file. Two addtional columns will be
                        appended to the orignal BED file with the last column
                        indicating "genes whose extended regulatory domain are
                        overlapped with the CpG", the 2nd last column
                        indicating "genes whose basal regulatory domain are
                        overlapped with the CpG". [required]
```	                        

### Input and output files

- BED3+ file specifying the C position. BED file must have at least 3 columns (chrom, start, end) 
```text
$ head test_01.bed6
#Chrom	Start	End	Name	Beta	Strand
chr1	10847	10848	cg26928153	0.8965	+
chr1	10849	10850	cg16269199	0.7915	+
chr1	15864	15865	cg13869341	0.9325	+
chr1	534241	534242	cg24669183	0.7941	+
chr1	564500	564501	cg26679879	0.3746	+
chr1	564503	564504	cg22519184	0.395	+
chr1	710096	710097	cg15560884	0.8106	+
chr1	714176	714177	cg01014490	0.0275	+
chr1	714620	714621	cg24063007	0.0368	+
```

- Reference gene model in standard BED-12 format
```text
$ head hg19.RefSeq.union.bed
chr1	11873	14409	DDX11L1	0	+	14409	14409	255,0,0	3	354,109,1189	0,739,1347
chr1	14361	29370	WASH7P	0	-	29370	29370	255,0,0	11	468,69,152,159,198,136,137,147,99,154,50	0,608,1434,2245,2496,2871,3244,3553,3906,10376,14959
chr1	17368	17436	MIR6859-1	0	-	17436	17436	255,0,0	1	68	0
chr1	17368	17436	MIR6859-2	0	-	17436	17436	255,0,0	1	68	0
chr1	34610	36081	FAM138A	0	-	36081	36081	255,0,0	3	564,205,361	0,666,1110
chr1	34610	36081	FAM138F	0	-	36081	36081	255,0,0	3	564,205,361	0,666,1110
chr1	69090	70008	OR4F5	0	+	69090	70008	255,0,0	1	918	0
chr1	134772	140566	LOC729737	0	-	140566	140566	255,0,0	3	4924,58,492	0,5017,5302
chr1	323891	328581	LOC100132062	0	+	328581	328581	255,0,0	3	169,58,4143	0,396,547
chr1	323891	328581	LOC100132287	0	+	328581	328581	255,0,0	3	169,58,4143	0,396,547
```

