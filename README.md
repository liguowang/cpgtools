## CpGtools -- Tools to analyze and visualize DNA methylation data

## Table of contents

- [Installation](#p1)
	- [Prerequisites](#p1.1)
	- [Python Dependencies](#p1.2)
	- [Install pip3](#p1.3)
	- [Install gamlss](#p1.4)
	- [Install CpGtools](#p1.5)
- [File and data format](#p2)
	- [BED format](#p2.1)
	- [Proportion value](#p2.2)
	- [The Beta-value](#p2.3)
	- [The M-value](#p2.4)
	- [Conversion between Beta- and M-value](#p2.5)
- [Usage information](#p3)

    |Program                              |Description                              
    |-------------------------------------|---------------------------------------------------------------------------------------------
    |[annotate_CpG.py](#p3.1)             |*Annotate the function of CpG by assigning it to gene's **basal** and **extended** regulatory domain*
    |[beta_m_conversion.py](#p3.2)        |*Convert Beta-value into M-value or vice versa*    
    |[beta_profile.py](#p3.3)             |*Calculate the average methylation level over genomic regions defined by genes (eg. exons, introns, intergenic regions)*
    |[chrom_distribution.py](#p3.4)       |*Calculates the distribution of CpG frequencies over chromosomes*
    |[dmc_bb.py](#p3.5)               	  |*Differential CpG analysis using **beta binomial regression** on proportion values* 
    |[dmc_fisher.py](#p3.6)               |*Differential CpG analysis using **Fisher's exact test** on proportion values* 
    |[dmc_glm.py](#p3.7)                  |*Differential CpG analysis using **linear model** on beta- or M-values* 
    |[dmc_logit.py](#p3.8)                |*Differential CpG analysis using **logistic regression model** on proportion values*
    |[dmc_nonparametric.py](#p3.9)        |*Differential CpG analysis using **MannWhitney U test** (2 groups comparison) or **KruskalWallis H-test** (3+ groups comparison) on beta- or M-values*
    |[dmc_ttest.py](#p3.10)                |*Differential CpG analysis using **T test** (2 groups comparison) or **ANOVA** (3+ groups comparison) on beta- or M-values*
    |[genomic_distribution_1.py](#p3.11)  |*Calculates the distribution of CpG frequencies over genomic regions defined by gene model* 
    |[genomic_distribution_2.py](#p3.12)  |*Calculates the distribution of CpG frequencies over genomic regions defined by user*
    |[methyl_logo.py](#p3.13)             |*Generate motif logo and motif matrices around cytosine*
    |[region_profile.py](#p3.14)          |*Calculate average methylation level for user specified genomic regions*
    |[region_stat.py](#p3.15)             |*Calculate basic statistics of CpGs located in each genomic region*
    |[trichotmize.py](#p3.16)             |*Trichotmize beta values into "methyl", "semimethyl" and "unmethyl" status using Gaussian Mixture Model* 

- [Comparison of differential CpG analysis tools](#p4)

- [Contact Information](#p5)

## <a name="p1"></a>Part 1: Installation


### <a name="p1.1"></a>Part 1.1: Prerequisites
CpGtools are written in [Python](https://www.python.org/). **Python3 (v3.5.x)**
is required to run all programs in CpGtools. Some programs also need **[R](https://www.r-project.org/)**
to generate graphs and fit linear and beta binomial models.  

- [Python 3](https://www.python.org/downloads/)
- [pip3](https://pip.pypa.io/en/stable/installing/)
- [R](https://www.r-project.org/)
- [gamlss](https://CRAN.R-project.org/package=gamlss) (Generalised Additive Models for Location Scale and Shape) R package

### <a name="p1.2"></a>Part 1.2: Python Dependencies
Note: You do NOT need to install these packages manually, as they will be automatically
installed when you use [pip3](https://pip.pypa.io/en/stable/installing/) to install CpGtools.

- [numpy](http://www.numpy.org/)
- [scipy](https://www.scipy.org/)
- [sklearn](https://www.scilearn.com/)
- [weblogo](https://pypi.org/project/weblogo/)
- [bx-python](https://github.com/bxlab/bx-python)

### <a name="p1.3"></a>Part 1.3: Install [pip3](https://pip.pypa.io/en/stable/installing/) (Skip this step if you already have it)

	$ curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
	$ python get-pip.py
	
Run the following code to check:

	$ pip3 --version
	pip 18.1 from /Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/pip (python 3.6)
	
	$ which pip3
	/Library/Frameworks/Python.framework/Versions/3.6/bin/pip
	
	$ which pip	
	/Library/Frameworks/Python.framework/Versions/3.6/bin/pip

Note that, in this case, **pip** is actually a soft link to **pip3**.

### <a name="p1.4"></a>Part 1.4: Install [gamlss](https://CRAN.R-project.org/package=gamlss)

If you don't have **[R](https://www.r-project.org/)**, please follow these [instructions](https://cran.r-project.org/doc/manuals/r-release/R-admin.html) to install.

Option-1: Install "gamlss" package from the R command line:

	> install.packages("gamlss")
	
Option-2: Install "gamlss" (use version 5.1-2 for example) from the Shell command line:
	
	$ wget https://cran.r-project.org/src/contrib/gamlss_5.1-2.tar.gz
	$ R CMD INSTALL gamlss_5.1-2.tar.gz

### <a name="p1.5"></a>Part 1.5: Install or Upgrade CpGtools
	
	$ pip3 install cpgtools		#install CpGtools and all its dependencies	
	$ pip3 install cpgtools --upgrade		#**upgrade** CpGtools and all its dependencies


## <a name="p2"></a>Part 2: File and data format 

### <a name="p2.1"></a>Part 2.1: BED format

BED (Browser Extensible Data) format is commonly used to describe *blocks of genome*.
The BED format consists of one line per feature, each containing 3-12 columns of data.
It is 0-based (meaning the first base of a chromosome is numbered 0). It is s left-open,
right-closed. For example, the bed entry **"chr1   10   15"** contains the 11-th, 12-th,
13-th, 14-th and 15-th bases of chromosome-1.

- **BED12** file (i.e. the standard BED file) which has 12 fields. Each row in this file describes a gene or an array of disconnected genomic regions. Details are described [here](https://genome.ucsc.edu/FAQ/FAQformat.html#format1). 
- **BED3** file only has the first three required fields (chrom, chromStart, chromEnd). Each row is used to represent a single genomic region where "score" and "strand" are not important. 
- **BED3+** file has at least three columns (chrom, chromStart, chromEnd). It could have additional columns, but these additional columns will be ignored.
- **BED6** file has the first six fields (chrom, chromStart, chromEnd, name, score, strand). Each row is used to represent a single genomic region and their associated scores, or in cases where "strand" information is important.  
- **BED6+** file has at least six columns (chrom, chromStart, chromEnd, name, score, stand). It could have additional columns, but these additional columns will be ignored.


### <a name="p2.2"></a>Part 2.2: proportion value
In [bisulfite sequencing](https://en.wikipedia.org/wiki/Bisulfite_sequencing) ([RRBS](https://en.wikipedia.org/wiki/Reduced_representation_bisulfite_sequencing) or [WGBS](https://en.wikipedia.org/wiki/Whole_genome_bisulfite_sequencing)), the methylation level of a particular CpG or region can be represented by a "proportion" vlaue. 
We define the proportion value as **two non-negative integers separated by comma (",")** with the first integer (*m*,  0 <= *m* <= *n*) representing 
"number of methylated reads" and the second integer (*n*, *n* >= 0) representing "number of total reads". for example:

```text
0,10 1,27 2,159		#Three proportions values indicated hypo-methylated loci 
7,7 17,19 30,34		#Three proportions values indicated hyper-methylated loci
```

### <a name="p2.3"></a>Part 2.3: Beta-value
The Beta-value is a value between 0 and 1, which can be interpreted as the approximation of the **percentage of methylation** for a given CpG or locus. 
One can convert *proportion value* into *beta value*, but not *vice versa*. In equation below, C is the "probe intensity" or "read count" of *methylated* allele, while U 
is the "probe intensity" or "read count" of *unmethylated* allele. 

![beta.png](https://github.com/liguowang/cpgtools/blob/master/img/beta.png)

### <a name="p2.4"></a>Part 2.4: M-value
The M-value is calculated as the log2 ratio of the probe intensities (or read counts) of *methylated* allele versus *unmethylated* allele.
In equation below, C is the "probe intensity" or "read count" of *methylated* allele, while U 
is the "probe intensity" or "read count" of *unmethylated* allele. w is the offset or pseudo count added to both denominator and numerator to avoid unexpected big changes and 
performing log transformation on zeros.

![M.png](https://github.com/liguowang/cpgtools/blob/master/img/M.png)

### <a name="p2.5"></a>Part 2.5: Conversion between Beta-value and M-value
The relationship between Beta-value and M-value is shown as equation and figure:

![beta_vs_M_curve.png](https://github.com/liguowang/cpgtools/blob/master/img/beta_vs_M_curve.png)


## <a name="p3"></a>Part 3: Usage Information


<a name="p3.1"></a>annotate_CpG.py
---
This program annotates CpGs by assigning them to their target genes. Follows the
"[Basal plus extension rules](http://great.stanford.edu/public/html/index.php)" used by [GREAT](http://great.stanford.edu/public/html/):

**Basal regulatory domain** is a user-defined genomic region around the TSS (transcription start site). By default,
from TSS upstream 5 Kb to TSS downstream 1 Kb is considered as the gene's *basal regulatory
domain*. When defining a gene's *basal regulatory domain*, the other nearby genes are
ignored (which means different genes' *basal regulatory domain* can be overlapped.)

**Extended regulatory domain** is a genomic region that is further extended from *basal regulatory domain* in both directions to the nearest gene's
basal regulatory domain but no more than the maximum extension (specified by '-e', default = 1000 kb) in one direction.	In other words, the "extension"
stops when it reaches other genes' "basal regulatory domain" or the extension limit, whichever comes first. 

*Basal regulatory domain* and *Extended regulatory domain* are illustrated in below diagram

![basal & extended regulatory domain](https://github.com/liguowang/cpgtools/blob/master/img/gene_domain.png)

#### Basic usage

```text

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        BED3+ file specifying the C position. BED3+ file could
                        be a regular text file or compressed file (*.gz,
                        *.bz2) or accessible url. [required]
  -r GENE_FILE, --refgene=GENE_FILE
                        Reference gene model in BED12 format
                        (https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
                        "One gene one transcript" is recommended. Since most
                        genes have multiple transcripts, one can collapse
                        multiple transcripts of the same gene into a single
                        super transcript or select the canonical transcript.
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
                        Prefix of the output file. Two additional columns will be
                        appended to the original BED file with the last column
                        indicating "genes whose extended regulatory domain are
                        overlapped with the CpG", the 2nd last column
                        indicating "genes whose basal regulatory domain are
                        overlapped with the CpG". [required]
```	                        

#### Input files
- BED3+ file specifying the C position. Download the test BED3 file [test_01.bed6](https://github.com/liguowang/cpgtools/blob/master/test/test_01.bed6)
- Reference gene model in BED12 format. Download the test BED12 file [hg19.RefSeq.union.bed](https://github.com/liguowang/cpgtools/blob/master/refgene/hg19.RefSeq.union.bed.gz)
- Human and mouse RefSeq gene bed files (with multiple transcripts of the same genes collapsed into a single super transcript) are available from [here](https://github.com/liguowang/cpgtools/blob/master/refgene/).

#### Output file
Two additional columns will be appended to the original BED file (-i):
- the last column contains genes whose **extended regulatory domain** are overlapped with the CpG
- the 2nd last column contains genes whose **basal regulatory domain** are overlapped with the CpG
- "//" indicates no genes are found

#### Example

```

$ python3 ../bin/annotate_CpG.py -r ../refgene/hg19.RefSeq.union.bed.gz  -i test_01.bed6 -o OUT1

@ 2018-12-07 12:49:21: Calculate basal regulatory domain from: "hg19.RefSeq.union.bed" ...
@ 2018-12-07 12:49:21: Calculate extended regulatory domain from: "hg19.RefSeq.union.bed" ...
@ 2018-12-07 12:49:22: Assigning CpG to gene ...

$ head OUT1.associated_genes.txt

#Chrom	Start	End	Name	Beta	Strand
chr1	10847	10848	cg26928153	0.8965	+	DDX11L1	//
chr1	10849	10850	cg16269199	0.7915	+	DDX11L1	//
chr1	15864	15865	cg13869341	0.9325	+	//	MIR6859-1;MIR6859-2
chr1	534241	534242	cg24669183	0.7941	+	//	OR4F29;OR4F3;LOC101928626;OR4F16
chr1	564500	564501	cg26679879	0.3746	+	LOC101928626	//
chr1	564503	564504	cg22519184	0.395	+	LOC101928626	//
chr1	710096	710097	cg15560884	0.8106	+	//	LOC100133331;LOC100288069
chr1	714176	714177	cg01014490	0.0275	+	LOC100288069	//
chr1	714620	714621	cg24063007	0.0368	+	LOC100288069	//

```


<a name="p3.2"></a>beta_m_conversion.py
---	
Convert Beta-value into M-value or vice versa.

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Data file with the 1st row containing sample IDs (must
                        be unique) and the 1st column containing CpG positions
                        or probe IDs (must be unique). This file can be a
                        regular text file or compressed file (*.gz, *.bz2) or
                        accessible url.
  -d DATA_TYPE, --dtype=DATA_TYPE
                        Data type either "Beta" or "M".
  -o OUT_FILE, --output=OUT_FILE
                        Output file.
```
#### Example

```
$ python3 ../bin/beta_m_conversion.py -i test_05_TwoGroup.tsv -d Beta -o test_05_Mvalue.txt

@ 2019-04-19 10:47:02: Convert Beta-value file "test_05_TwoGroup.tsv" into M-value file "test_05_Mvalue.txt" ...

```





<a name="p3.3"></a>beta_profile.py
---	
beta_profile.py calculates the average methylation level (i.e. average beta value) across
gene body (including: 5'UTR exon, CDS exon, 3'UTR exon, first intron, internal intron, last
intron,  up-stream intergenic and down-stream intergenic regions).

#### Basic usage

```text

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        BED6+ file specifying the C position. This BED file
                        should have at least 6 columns (Chrom, ChromStart,
                        ChromeEnd, Name, Beta_value, Strand). BED6+ file can
                        be a regular text file or compressed file (*.gz,
                        *.bz2) or accessible url.
  -r GENE_FILE, --refgene=GENE_FILE
                        Reference gene model in standard BED12 format
                        (https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
                        "Strand" column must exist in order to decide 5' and
                        3' UTRs, up- and down-stream intergenic regions.
  -d DOWNSTREAM_SIZE, --downstream=DOWNSTREAM_SIZE
                        Size of down-stream genomic region added to gene.
                        default=2000 (bp)
  -u UPSTREAM_SIZE, --upstream=UPSTREAM_SIZE
                        Size of up-stream genomic region added to gene.
                        default=2000 (bp)
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.                        
```

#### Input files
- BED6+ file specifying the C position. Download test file [test_02.bed6.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_02.bed6.gz)
- Reference gene model in BED12 format. Download test file [hg19.RefSeq.union.bed](https://github.com/liguowang/cpgtools/blob/master/refgene/hg19.RefSeq.union.bed.gz)

#### Example

```
$ python3 ../bin/beta_profile.py -r ../refgene/hg19.RefSeq.union.bed.gz -i test_02.bed6 -o OUT2

@ 2018-12-07 13:43:00: Reading CpG file: "test_02.bed6.gz"
@ 2018-12-07 13:43:09: Reading reference gene model: "hg19.RefSeq.union.bed"
@ 2018-12-07 13:43:09: Process upstream regions ...
@ 2018-12-07 13:43:10: Process 5' UTR exons ...
@ 2018-12-07 13:43:10: Process Coding exons ...
@ 2018-12-07 13:43:11: Process first introns ...
@ 2018-12-07 13:43:12: Process internal introns ...
@ 2018-12-07 13:43:13: Process last introns ...
@ 2018-12-07 13:43:14: Process 3' UTR exons ...
@ 2018-12-07 13:43:15: Process downstream regions ...

```

#### Output
- The red curve represents average profile of beta values, aggregated from gene regions defined in BED file (-r). 
- Upstream: intergenic region (defined by '-u' ) before TSS (transcription start site)
- Downstream: intergenic region (defined by '-d') after TES (transcription end site)
![beta_profile.png](https://github.com/liguowang/cpgtools/blob/master/img/beta_profile.png)


<a name="p3.4"></a>chrom_distribution.py
---	
This program calculates the distribution of CpG over chromosomes.

#### Basic usage

```text

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILES, --input-files=INPUT_FILES
                        Input CpG file(s) in BED3+ format. Multiple BED files
                        should be separated by "," (eg: "-i
                        file_1.bed,file_2.bed,file_3.bed"). BED file can be
                        regular or compressed by 'gzip' or 'bz'. The barplot
                        figures will NOT be generated if you provide more than
                        12 samples (bed files). [required]
  -n FILE_NAMES, --names=FILE_NAMES
                        Shorter and meaningful names to label samples. Should
                        be separated by "," and match CpG BED files in number.
                        If not provided, basenames of CpG BED files will be
                        used to label samples. [optional]
  -s CHROM_SIZE, --chrom-size=CHROM_SIZE
                        Chromosome size file. Tab or space separated text file
                        with two columns: the first column is chromosome name/ID,
                        the second column is chromosome size. This file will
                        determine: (1) which chromosomes are included in the
                        final barplots, so do NOT include 'unplaced',
                        'alternative' contigs in this file. (2) The order of
                        chromosomes in the final barplots.  [required]
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file. [required]

```
#### Input files

- BED3+ files: [test_03a.bed3.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_03a.bed3.gz), [test_03b.bed3.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_03b.bed3.gz)
- chromosome size file: [hg19.chrom.sizes](https://github.com/liguowang/cpgtools/blob/master/test/hg19.chrom.sizes)
#### Example

```text

$ python3 ../bin/chrom_distribution.py -i test_03a.bed3.gz,test_03b.bed3.gz -n 450K,850K -s hg19.chrom.sizes -o chromDist

```

#### Output files

1. Total CpG count per chromsome 
![chromDist.CpG_total.png](https://github.com/liguowang/cpgtools/blob/master/img/chromDist.CpG_total.png) 
2. CpG percent on each chromosome (normalized to total CpGs)
![chromDist.CpG_percent.png](https://github.com/liguowang/cpgtools/blob/master/img/chromDist.CpG_percent.png)
3. CpG per Mb (normalized to chromsome size)
![chromDist.CpG_perMb.png](https://github.com/liguowang/cpgtools/blob/master/img/chromDist.CpG_perMb.png)


<a name="p3.5"></a>dmc_bb.py
---
This program performs differential CpG analysis using **"beta binomial (BB)"** or **"zero inflated beta binomial model (ZIBB)"** on [proportion values](#p2.2). It allows for covariable analysis.

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Data file containing methylation proportions
                        (represented by "methyl_count,total_count", eg.
                        "20,30") with the 1st row containing sample IDs (must
                        be unique) and the 1st column containing CpG positions
                        or probe IDs (must be unique). This file can be a
                        regular text file or compressed file (*.gz, *.bz2) or
                        accessible url.
  -g GROUP_FILE, --group=GROUP_FILE
                        Group file defining the biological groups of each
                        sample as well as other covariables such as gender,
                        age.  Sample IDs should match to the "Data file".
  -f FAMILY_FUNC, --family=FAMILY_FUNC
                        A gamlss (https://cran.r-project.org/web/packages/gaml
                        ss/index.html) family object. Can be integer 1 or 2
                        with 1 = "BB (beta binomial)", 2 = "ZIBB (zero
                        inflated beta binomial)" or 3 = "ZABB (zero adjusted
                        beta binomial)". Default=1.
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
```                        

#### Example
```text
$ python3 ../bin/dmc_bb.py -i test_04_TwoGroup.tsv.gz -g test_04_TwoGroup.grp.csv.gz -o OUT_4
$ python3 ../bin/dmc_bb.py -i test_04_TwoGroup.tsv.gz -g test_04_TwoGroup.grp.csv.gz -o OUT_4
```

#### Output file
Additional columns (pvalue and coefficient) will be appended to the original data file. In the example above, 
4 additional columns were added to "test_04_TwoGroup.tsv":

- survival.pval
- Sex.pval
- survival.coef
- Sex.coef

<a name="p3.6"></a>dmc_fisher.py
---
This program performs differential CpG analysis using **Fisher exact test** on [proportion value](#p2.2).

 * applies to two sample comparison with no biological/technical replicates
 * if biological/technical replicates are provided, *methyl reads* and *total reads* of all replicates will be merged (i.e. ignores biological/technical variations)

Input file format:

```text
cgID        sample_1    sample_2	...
CpG_1       129,170     166,178	...
CpG_2       24,77       67,99	...
...
```

number before "," indicates *number of methyl reads*
number after "," indicates *number of total reads*

3 columns ("Odds ratio", "pvalue" and "FDR adjusted pvalue") will append to this table.

- pvalue is two-tailed (same for other methods)
- We used [Benjamini-Hochberg Procedure](https://www.statisticshowto.datasciencecentral.com/benjamini-hochberg-procedure/) for multiple test correction (same for other methods)

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Data file containing methylation proportions
                        (represented by "methyl_count,total_count", eg.
                        "20,30") with the 1st row containing sample IDs (must
                        be unique) and the 1st column containing CpG positions
                        or probe IDs (must be unique). This file can be a
                        regular text file or compressed file (*.gz, *.bz2) or
                        accessible url.
  -g GROUP_FILE, --group=GROUP_FILE
                        Group file defining the biological group of each
                        sample. It is a comma-separated two columns file with
                        the 1st column containing sample IDs, and the 2nd
                        column containing group IDs.  It must have a header
                        row. Sample IDs should match to the "Data file".
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
```

#### Input files

- Data file:[test_09.tsv](https://github.com/liguowang/cpgtools/blob/master/test/test_09.tsv)
- Group file:[test_09.grp.csv](https://github.com/liguowang/cpgtools/blob/master/test/test_09.grp.csv)


#### Example

```text
$ python3 ../bin/dmc_fisher.py  -i test_09.tsv -g test_09.grp.csv -o test_09

@ 2018-12-13 14:40:28: Read group file "test_09.grp.csv" ...
	Group 1 has 2 samples:
		LTS_MCR-1008,LTS_MCR-1035
	Group 2 has 2 samples:
		STS_MCR-1021,STS_MCR-1251
@ 2018-12-13 14:40:48: Perfrom Benjamini-Hochberg (aka FDR) correction ...
@ 2018-12-13 14:40:48: Writing to test_09.pval.txt


$ head test_09.pval.txt

ID	LTS_MCR-1008	LTS_MCR-1035	STS_MCR-1021	STS_MCR-1251	OddsRatio	pval	adj.pval
chr10:100011340	12,14	26,37	0,18	10,24	9.353846153846154	1.2116597355208375e-06	6.343768248800197e-05
chr10:100011341	0,21	0,54	0,26	0,19	nan	1.0	1.0
chr10:100011387	0,14	0,40	0,20	0,24	nan	1.0	1.0
chr10:100011388	18,18	47,54	19,23	18,19	1.2548262548262548	0.7574366471769988	1.0
chr10:100026933	16,30	28,55	7,40	13,19	2.0926829268292684	0.04119183894184185	0.2617016451197068
chr10:100026991	0,30	0,48	0,40	0,19	nan	1.0	1.0
chr10:100027909	2,77	1,66	2,64	0,46	1.1571428571428573	1.0	1.0
chr10:100027910	0,34	0,49	0,49	0,40	nan	1.0	1.0
chr10:100027919	0,76	0,66	2,58	0,44	0.0	0.17375025298519042	0.6757824934416998
		
```       

<a name="p3.7"></a>dmc_glm.py
---
This program performs differential CpG analysis using **[generalized liner model](https://en.wikipedia.org/wiki/Generalized_linear_model)** based on
[beta values](#p2.3). It allows for covariable analysis.

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Data file containing beta values with the 1st row
                        containing sample IDs (must be unique) and the 1st
                        column containing CpG positions or probe IDs (must be
                        unique). This file can be regular or compressed by
                        'gzip' or 'bz'.
  -g GROUP_FILE, --group=GROUP_FILE
                        Group file defining the biological group of each
                        sample as well as other covariables such as gender,
                        age.  Sample IDs should match to the "Data file".
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.

```

##### Input files

- Data file: [test_10_TwoGroup.tsv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_10_TwoGroup.tsv.gz)
- Group file: [test_10_TwoGroup.grp.csv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_10_TwoGroup.grp.csv.gz)

#### Output files
Additional columns (pvalue and coefficient) will be appended to the original data file.

<a name="p3.8"></a>dmc_logit.py
---
This program performs differential CpG analysis using [logistic regression](https://en.wikipedia.org/wiki/Logistic_regression) model based on
[proportion values](#p2.2). It allows for covariable analysis.
Users can choose to use "binomial" or "quasibinomial" to model the data. According to "glm" documentation, "The *quasibinomial* family differs from the *binomial* family only in that the dispersion parameter is not fixed at one, so it can model over-dispersion".

#### Basic usage

```text

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Data file containing methylation proportions
                        (represented by "methyl_count,total_count", eg.
                        "20,30") with the 1st row containing sample IDs (must
                        be unique) and the 1st column containing CpG positions
                        or probe IDs (must be unique). This file can be a
                        regular text file or compressed file (*.gz, *.bz2) or
                        accessible url..
  -g GROUP_FILE, --group=GROUP_FILE
                        Group file defining the biological group of each
                        sample as well as other covariables such as gender,
                        age.  Sample IDs should match to the "Data file".
  -f FAMILY_FUNC, --family=FAMILY_FUNC
                        Error distribution and link function to be used in the
                        GLM model. Can be integer 1 or 2 with 1 = "binomial"
                        and 2 = "quasibinomial". Default=1.
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
```

#### Input files

- Data file: [test_04_TwoGroup.tsv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_04_TwoGroup.tsv.gz)
- Group file: [test_04_TwoGroup.grp.csv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_04_TwoGroup.grp.csv.gz)

#### Example
```text
$ python3 ../bin/dmc_logit.py -i test_04_TwoGroup.tsv.gz -g test_04_TwoGroup.grp.csv.gz -o OUT_4
```

#### Output file
Additional columns (pvalue and coefficient) will be appended to the original data file. In the example above, 
4 additional columns were added to "test_04_TwoGroup.tsv":

- survival.pval
- Sex.pval
- survival.coef
- Sex.coef

<a name="p3.9"></a>dmc_nonparametric.py
---
This program performs differential CpG analysis based on [beta values](#p2.3).
- use [Mann-Whitney U test](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html) for two group comparison.
- use [Kruskal-Wallis H-test](https://en.wikipedia.org/wiki/Kruskal%E2%80%93Wallis_one-way_analysis_of_variance) for multiple groups comparison.

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Data file containing beta values with the 1st row
                        containing sample IDs (must be unique) and the 1st
                        column containing CpG positions or probe IDs (must be
                        unique). Except for the 1st row and 1st column, any
                        non-numerical values will be considered as "missing
                        values" and ignored. This file can be regular or
                        compressed by 'gzip' or 'bz'.
  -g GROUP_FILE, --group=GROUP_FILE
                        Group file defining the biological group of each
                        samples. It is a comma-separated two columns file with
                        the 1st column containing sample IDs, and the 2nd
                        column containing group IDs. It must have a header
                        row. Sample IDs should match to the "Data file". Note:
                        automatically switch to use  Kruskal-Wallis H-test if
                        more than 2 groups were defined in this file.
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
```

#### Input files
- [test_05_TwoGroup.tsv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_05_TwoGroup.tsv.gz) 
- [test_05_TwoGroup.grp.csv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_05_TwoGroup.grp.csv.gz)      
- [test_06_ThreeGroup.tsv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_06_ThreeGroup.tsv.gz)    
- [test_06_ThreeGroup.grp.csv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_06_ThreeGroup.grp.csv.gz)          

#### Example
```text
$ python3 ../bin/dmc_nonparametric.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv.gz -o OUT_05
@ 2018-12-11 11:17:42: Read group file "test_05_TwoGroup.grp.csv.gz" ...
	Group 1 has 10 samples:
		Normal_01,Normal_02,Normal_03,Normal_04,Normal_05,Normal_06,Normal_07,Normal_08,Normal_09,Normal_10
	Group 2 has 10 samples:
		CirrHCV_01,CirrHCV_02,CirrHCV_03,CirrHCV_04,CirrHCV_05,CirrHCV_06,CirrHCV_07,CirrHCV_08,CirrHCV_09,CirrHCV_10
@ 2018-12-11 11:17:42: Perfrom Mann-Whitney rank test of two samples ...
@ 2018-12-11 11:17:45: Perfrom Benjamini-Hochberg (aka FDR) correction ...
@ 2018-12-11 11:17:46: Writing to OUT_05.pval.txt



$ python3 ../bin/dmc_nonparametric.py -i test_06_ThreeGroup.tsv.gz -g test_06_ThreeGroup.grp.csv.gz -o OUT_06
@ 2018-12-11 11:18:34: Read group file "test_06_ThreeGroup.grp.csv.gz" ...
	Group 1 has 10 samples:
		Normal_01,Normal_02,Normal_03,Normal_04,Normal_05,Normal_06,Normal_07,Normal_08,Normal_09,Normal_10
	Group 2 has 10 samples:
		CirrHCV_01,CirrHCV_02,CirrHCV_03,CirrHCV_04,CirrHCV_05,CirrHCV_06,CirrHCV_07,CirrHCV_08,CirrHCV_09,CirrHCV_10
	Group 3 has 10 samples:
		HCCHCV_01,HCCHCV_02,HCCHCV_03,HCCHCV_04,HCCHCV_05,HCCHCV_06,HCCHCV_07,HCCHCV_08,HCCHCV_09,HCCHCV_10
@ 2018-12-11 11:18:34: Perfrom Kruskal-Wallis H-test ...
@ 2018-12-11 11:18:40: Perfrom Benjamini-Hochberg (aka FDR) correction ...
@ 2018-12-11 11:18:40: Writing to OUT_06.pval.txt

```   
#### Output file
Additional two columns ("pval", and "adj.pval") will be appended to the orignal data file.


<a name="p3.10"></a>dmc_ttest.py
---
This program performs differential CpG analysis based on [beta values](#p2.3).

* uses [Student's t-test](https://en.wikipedia.org/wiki/Student%27s_t-test) for two group comparison.
* uses [ANOVA](https://en.wikipedia.org/wiki/Analysis_of_variance) for multiple groups comparison.


#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Data file containing beta values with the 1st row
                        containing sample IDs (must be unique) and the 1st
                        column containing CpG positions or probe IDs (must be
                        unique). Except for the 1st row and 1st column, any
                        non-numerical values will be considered as "missing
                        values" and ignored. This file can be regular or
                        compressed by 'gzip' or 'bz'.
  -g GROUP_FILE, --group=GROUP_FILE
                        Group file defining the biological groups of each
                        sample. It is a comma-separated 2 columns file with
                        the 1st column containing sample IDs, and the 2nd
                        column containing group IDs.  It must have a header
                        row. Sample IDs shoud match to the "Data file". Note:
                        automatically switch to use ANOVA if more than 2
                        groups were defined in this file.
  -p, --paired          If '-p/--paired' flag was specified, use paired t-test
                        which requires the equal number of samples in both
                        groups. Paired sampels are matched by the order. This
                        option will be ignored for multiple group analysis.
  -w, --welch           If '-w/--welch' flag was specified, using Welch's
                        t-test which does not assume the two samples have
                        equal variance.  If omitted, use standard two-sample
                        t-test (i.e. assuming the two samples have equal
                        variance). This option will be ignored for paired
                        t-test and multiple group analysis.
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.    
```

#### Example
```text
$ python3 ../bin/dmc_ttest.py -i test_05_TwoGroup.tsv.gz -g test_05_TwoGroup.grp.csv.gz -o OUT_05

@ 2018-12-11 12:36:48: Read group file "test_05_TwoGroup.grp.csv.gz" ...
	Group 1 has 10 samples:
		Normal_01,Normal_02,Normal_03,Normal_04,Normal_05,Normal_06,Normal_07,Normal_08,Normal_09,Normal_10
	Group 2 has 10 samples:
		CirrHCV_01,CirrHCV_02,CirrHCV_03,CirrHCV_04,CirrHCV_05,CirrHCV_06,CirrHCV_07,CirrHCV_08,CirrHCV_09,CirrHCV_10
@ 2018-12-11 12:36:48: Perfrom standard t-test of two independent samples ...

@ 2018-12-11 12:36:52: Perfrom Benjamini-Hochberg (aka FDR) correction ...
@ 2018-12-11 12:36:52: Writing to OUT_05.pval.txt



$ python3 ../bin/dmc_ttest.py -i test_06_ThreeGroup.tsv.gz -g test_06_ThreeGroup.grp.csv.gz -o OUT_06

@ 2018-12-11 12:37:43: Read group file "test_06_ThreeGroup.grp.csv.gz" ...
	Group 1 has 10 samples:
		Normal_01,Normal_02,Normal_03,Normal_04,Normal_05,Normal_06,Normal_07,Normal_08,Normal_09,Normal_10
	Group 2 has 10 samples:
		CirrHCV_01,CirrHCV_02,CirrHCV_03,CirrHCV_04,CirrHCV_05,CirrHCV_06,CirrHCV_07,CirrHCV_08,CirrHCV_09,CirrHCV_10
	Group 3 has 10 samples:
		HCCHCV_01,HCCHCV_02,HCCHCV_03,HCCHCV_04,HCCHCV_05,HCCHCV_06,HCCHCV_07,HCCHCV_08,HCCHCV_09,HCCHCV_10
@ 2018-12-11 12:37:43: Perfrom ANOVA ...
@ 2018-12-11 12:37:45: Perfrom Benjamini-Hochberg (aka FDR) correction ...
@ 2018-12-11 12:37:45: Writing to OUT_06.pval.txt
```
#### Output file
Additional two columns ("pval", and "adj.pval") will be appended to the orignal data file.


              

<a name="p3.11"></a>genomic_distribution_1.py
----
This program counts number of CpGs falling into genomic regions around **genes**. The genomic region around a particular gene can be divided into 5 groups:

1. Coding exons
2. UTR exons
3. Introns
4. Upstream intergenic regions (regions upstream of TSS)
5. Downsteam intergenic regions (regions downstream of TES)

Please note, a particular genomic region can be assigned to different groups defined above,
because most genes have multiple transcripts, and different genes could overlap on the
genome. For example, an exon of one gene could be located in (or partially overlapped with) an intron of another gene. To address
this ambiguity issue, we define the following priority order:

*Coding exons* > *UTR exons* > *Introns* > *Upstream intergenic regions* > *Downsteam intergenic regions*

Higher-priority group overrides the low-priority group. For example, if a certain part
of an **intron** is overlapped with **exon** of other transcripts/genes, the overlapped part will
be considered as exon (i.e. removed from intron) since "exon" has higher priority.

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        BED file specifying the methylated C position. This
                        BED file should have at least 3 columns (Chrom,
                        ChromStart, ChromeEnd).  Note: the first base in a
                        chromosome is numbered 0. BED file can be regular or
                        compressed by 'gzip' or 'bz'.
  -r GENE_FILE, --refgene=GENE_FILE
                        Reference gene model in standard BED-12 format
                        (https://genome.ucsc.edu/FAQ/FAQformat.html#format1).
  -d DOWNSTREAM_SIZE, --downstream=DOWNSTREAM_SIZE
                        Size of down-stream intergenic region w.r.t. TES
                        (transcription end site). default=2000 (bp)
  -u UPSTREAM_SIZE, --upstream=UPSTREAM_SIZE
                        Size of up-stream intergenic region w.r.t. TSS
                        (transcription start site). default=2000 (bp)
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.

```

#### Example
```text

$ python3 ../bin/genomic_distribution_1.py -i test_03b.bed3.gz -r hg19.RefSeq.union.bed -o OUT_7

```

#### Output files
```text
$ cat OUT_7.tsv

Priority_order	Name	Number_of_regions	Size_of_regions(bp)	CpG_raw_count	CpG_count_per_KB
0	Coding exons	204685	39119881	65488	1.674033722137345
1	UTR exons	69937	38385741	61510	1.6024179395156133
2	Introns	214085	1228745034	329012	0.26776262845103793
3	Upstream of TSS	20507	37014855	120353	3.251478359161477
4	Downstream of TES	18790	35709088	10999	0.3080168275370124

The barplot "OUT_7.pdf" was also generated.
```         
![Genomic distribution.png](https://github.com/liguowang/cpgtools/blob/master/img/genomic_dist1.png)  

<a name="p3.12"></a>genomic_distribution_2.py
----
This program counts number of CpGs falling into genomic regions defined by **users**.
A maximum of 10 BED files (defining 10 sets of genomic regions) can be analyzed.

Please note:
The order of BED files is important (i.e. considered as "priority order"). Overlapped
genomic regions will be retained only in the **set** with the highest priority and removed from
**all the other sets** that have lower priorities.  For example, users provided 3 BED files via
"-i promoters.bed,enhancers.bed,intergenic.bed", then if an enhancer region is overlapped
with promoters, **the overlapped part** will be removed from "enhancers.bed".

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i CPG_FILE, --cpg=CPG_FILE
                        BED file specifying the methylated C position. This
                        BED file should have at least 3 columns (Chrom,
                        ChromStart, ChromeEnd).  Note: the first base in a
                        chromosome is numbered 0. BED file can be regular or
                        compressed by 'gzip' or 'bz'.
  -b BED_FILES, --bed=BED_FILES
                        List of BED files specifying the genomic regions.
                        Note: (1) This program can only analyze a maximum of
                        10 BED files. (2) BED files should be separated by
                        comma (eg. " -i
                        promoters.bed,enhancers.bed,intergenic.bed"). (3) The
                        *order* of BED files is used to determine the
                        *priority* of BED files, and overlapped genomic
                        regions will be kept only in the BED file of the
                        highest priority and removed from BED files of lower
                        priority. For example, if an enhancer region is
                        overlapped with promoters, the *overlapped part* will
                        be removed from "enhancers.bed". (4) Each BED file
                        should have at least 3 columns (Chrom, ChromStart,
                        ChromeEnd), and the first base in a chromosome is
                        numbered 0. (5) BED files can be regular or compressed
                        by 'gzip' or 'bz'.
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
```

#### Input files
- [hg19_H3K4me3.bed4](https://github.com/liguowang/cpgtools/blob/master/test/hg19_H3K4me3.bed4)
- [hg19_CGI.bed4](https://github.com/liguowang/cpgtools/blob/master/test/hg19_CGI.bed4)
- [hg19_H3K27ac_with_H3K4me1.bed4](https://github.com/liguowang/cpgtools/blob/master/test/hg19_H3K27ac_with_H3K4me1.bed4)
- [hg19_H3K27me3.bed4](https://github.com/liguowang/cpgtools/blob/master/test/hg19_H3K27me3.bed4)
            
#### Example
```text
$ python3 ../bin/genomic_distribution_2.py -i test_03b.bed3.gz  -b  hg19_H3K4me3.bed4,hg19_CGI.bed4,hg19_H3K27ac_with_H3K4me1.bed4,hg19_H3K27me3.bed4 -o OUT_8

```

#### Output files

OUT_8.txt

|Priority_order  |Name                             |Number_of_regions       |Size_of_regions(bp)     |CpG_raw_count   |CpG_count_per_KB
|----------------|---------------------------------|------------------------|------------------------|----------------|-------------------
|0               |hg19_H3K4me3.bed4                |88439                   |215137961               |253264          |1.177216697707756
|1               |hg19_CGI.bed4                    |13637                   |6739408                 |42559           |6.314946357306161
|2               |hg19_H3K27ac_with_H3K4me1.bed4   |49579                   |83443517                |51436           |0.6164169710152557
|3               |hg19_H3K27me3.bed4               |137573                  |276362203               |60922           |0.22044259069681826      

OUT_8.pdf
![Genomic distribution2.png](https://github.com/liguowang/cpgtools/blob/master/img/genome_dist2.png)  

<a name="p3.13"></a>methyl_logo.py
----
This program generates DNA sequence logo around methylated Cs in 3 steps:

1. Extract genomic sequences around methylated C postion
2. Generate [motif matrices](https://en.wikipedia.org/wiki/Position_weight_matrix) include:
	- position frequency matrix (PFM)
	- position probability matrix (PPM)
	- position weight matrix (PWM)
	- [MEME](http://meme-suite.org/doc/meme-format.html) format matrix
	- [Jaspar](http://jaspar.genereg.net/) format matrix
3. Generate motif logo using [weblogo](https://github.com/WebLogo/weblogo)

Note: **input BED file** must has strand information. 

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        BED file of the methylated C position. This
                        BED file should have at least 6 columns (Chrom,
                        ChromStart, ChromeEnd, name, score, strand).  Note:
                        Must provide correct *strand* information. BED file
                        can be regular or compressed by 'gzip' or 'bz'.
  -r GENOME_FILE, --refgenome=GENOME_FILE
                        Reference genome seqeunces in FASTA format.
  -e EXTEND_SIZE, --extend=EXTEND_SIZE
                        Number of bases extended to up- and down-stream.
                        default=5 (bp)
  -n MOTIF_NAME, --name=MOTIF_NAME
                        Motif name. default=motif
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
```                        
#### Input file
[test_07_450K_CH.bed](https://github.com/liguowang/cpgtools/blob/master/test/test_07_450K_CH.bed)

#### Example
```text
$ python3 ../bin/methyl_logo.py -i test_07_450K_CH.bed -r /database/hg19.fa -o OUT_9
```

#### Output

Motif logo

![methyl_logo.png](https://github.com/liguowang/cpgtools/blob/master/img/methyl_logo.png)

Motif marices
```text
$ cat OUT_9.pfm
Base	1	2	3	4	5	6	7	8	9	10	11
A	969.0	979.0	299.0	211.0	2838.0	26.0	2891.0	262.0	81.0	1122.0	1083.0
C	412.0	447.0	831.0	190.0	22.0	2898.0	1.0	168.0	1291.0	1242.0	860.0
G	876.0	571.0	194.0	35.0	48.0	7.0	9.0	2477.0	739.0	58.0	462.0
T	675.0	935.0	1608.0	2496.0	24.0	1.0	31.0	25.0	821.0	510.0	527.0

$ cat OUT_9.meme
----------------------------------------
motif position-specific probability matrix
----------------------------------------
letter-probability matrix: alength= 4 w= 11 nsites= 2932
 0.3304691762138571 0.1405482815057283 0.2987588652482269 0.23022367703218768
 0.33387888707037644 0.1524822695035461 0.19476268412438624 0.3188761593016912
 0.10201854882705945 0.2834151663938898 0.0662165848336061 0.5483496999454446
 0.07201309328968904 0.06485270049099837 0.012002182214948174 0.8511320240043645
 0.967744135297327 0.007569558101472996 0.016434806328423354 0.00825150027277687
 0.008933442444080744 0.9882024004364431 0.0024549918166939452 0.00040916530278232413
 0.9858156028368796 0.00040916530278232413 0.0031369339879978183 0.010638297872340429
 0.0894026186579378 0.057351336606655756 0.8446535733769777 0.008592471358428807
 0.02768685215493726 0.4402618657937807 0.2520458265139116 0.28000545553737044
 0.3826377523186034 0.42355428259683586 0.01984451718494272 0.17396344789961812
 0.36933987997817785 0.29330332787779595 0.15759683578832515 0.17975995635570105
```

<a name="p3.14"></a>region_profile.py
----
This program calculates the overall methylation level (i.e. average beta value) over
particular genomic regions (eg. promoters, TF bindings). 

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        BED6 file of the C position. This BED file
                        should have at least 6 columns (Chrom, ChromStart,
                        ChromeEnd, Name, Beta_value, Strand).  Note: the first
                        base in a chromosome is numbered 0. BED file can be
                        regular or compressed by 'gzip' or 'bz'.
  -r REGION_FILE, --region=REGION_FILE
                        BED file of genomic regions. This BED file
                        should have at least three columns (Chrom, ChromStart,
                        ChromeEnd). If the 6-th column does not exist, all
                        regions will be considered as on "+" strand.
  -d DOWNSTREAM_SIZE, --downstream=DOWNSTREAM_SIZE
                        Size of extension to downstream. default=2000 (bp)
  -u UPSTREAM_SIZE, --upstream=UPSTREAM_SIZE
                        Size of extension to upstream. default=2000 (bp)
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of output file.
```                        
#### Example
```text
python3 ../bin/region_profile.py -i test_02.bed6.gz -r hg19.RefSeq.union.1Kpromoter.bed -o OUT_10
```

#### Output
- Average methylation profile over "upstream"(left), "user-specified" (middle), and "down-stream"(right) genomic regions. 
- Note: All groups ("upstream", "user-specified", and "down-stream") were scaled into 100 nt long, regardless of their actual genomic sizes.   

![region_profile.png](https://github.com/liguowang/cpgtools/blob/master/img/region_profile.png)

<a name="p3.15"></a>region_stat.py
----
This program gives basic statistics for each genomic region. Append 6 columns to the input BED file:

1. Number of CpGs detected in the genomic region
2. Min methylation level
3. Max methylation level
4. Average methylation level across all CpGs
5. Median methylation level across all CpGs
6. Standard deviation

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        BED file specifying the methylated C position. This
                        BED file should have at least 6 columns (Chrom,
                        ChromStart, ChromeEnd, Name, Beta_value, Strand).
                        Note: the first base in a chromosome is numbered 0.
                        BED file can be regular or compressed by 'gzip' or
                        'bz'.
  -r REGION_FILE, --region=REGION_FILE
                        BED file specificy genomic regions. This BED file
                        should have at least 3 columns (Chrom, ChromStart,
                        ChromeEnd).
  -o OUT_FILE, --output=OUT_FILE
                        Prefix of the output file.
```

#### Example
```text

$ python3 ../bin/region_stat.py -i test_02.bed6.gz -r hg19.RefSeq.union.1Kpromoter.bed -o OUT_11

$ grep -v "NA" OUT_11.txt |head
chr1	563888	564889	LOC101928626	0	-	2	0.3746	0.395	0.3848	0.3848	0.0144
chr1	752250	753251	FAM87B	0	+	2	0.4641	0.6485	0.5563	0.5563	0.1304
chr1	811681	812682	FAM41C	0	-	3	0.7784	0.8494	0.8228	0.8407	0.0387
chr1	854571	855572	LOC100130417	0	-	11	0.4197	0.9319	0.8018	0.8744	0.1602
chr1	860620	861621	SAMD11	0	+	8	0.0252	0.5	0.1595	0.0505	0.1972
chr1	894178	895179	NOC2L	0	-	4	0.0571	0.234	0.137	0.1285	0.0755
chr1	895466	896467	KLHL17	0	+	7	0.0339	0.4091	0.1589	0.0924	0.1405
chr1	901376	902377	PLEKHN1	0	+	7	0.0294	0.7146	0.2418	0.1511	0.2615
chr1	916996	917997	PERM1	0	-	8	0.2443	0.6882	0.4772	0.481	0.1363
chr1	935051	936052	HES4	0	-	3	0.0863	0.1507	0.1145	0.1066	0.0329


```     

<a name="p3.16"></a>trichotmize.py
----
This program uses Gaussian Mixture model (GMM) to trichotmize each CpG into 4 status:
 * Un-methylated (labeled as "0" in result file)
 * Semi-methylated (labeled as "1" in result file)
 * Full-methylated (labeled as "2" in result file)
 * Unknown

Basically, GMM will first calculate probability p (p = {p0 + p1 + p2}, p0 + p1 + p2 = 1) for each CpG based on its beta value
- p0: the probability that the CpG is un-methylated
- p1: the probability that the CpG is semi-methylated
- p2: the probability that the CpG is full-methylated

The classification will be made using rules: 

- if p0 == max(p): un-methylated
- if p2 == max(p): full-methylated
- if p1 == max(p):
	- if p1 <= prob_cutoff: semi-methylated
	- else: unknown
 
#### Input files
[test_08.tsv.gz](https://github.com/liguowang/cpgtools/blob/master/test/test_08.tsv.gz) 

#### Basic usage
```text
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -i INPUT_FILE, --input-file=INPUT_FILE
                        Input plain text file containing beta values with the
                        1st row containing sample IDs (must be unique) and the
                        1st column containing probe IDs (must be unique).
                        [required]
  -c PROB_CUTOFF, --prob-cut=PROB_CUTOFF
                        Probability cutoff to assign a probe into "semi-
                        methylated" class. default=0.9999
  -r, --report          Presense of this flag renders program to generate
                        "summary_report.txt" file.
  -s RANDOM_STATE, --seed=RANDOM_STATE
                        Random_state is the seed used by the random number
                        generator. You get exactly the same results when
                        running multiple times with the same random_state.
                        default=1
``` 
#### Example
```text
$ python3 ../bin/trichotmize.py -r -i test_08.tsv.gz

@ 2018-12-11 16:25:39: Reading test_08.tsv.gz ...
	Total samples: 5
	Total probes: 386958
@ 2018-12-11 16:25:42: Building Bayesian Gaussian Mixture model for subject: TCGA-BC-A10Q ...
@ 2018-12-11 16:25:50: Building Bayesian Gaussian Mixture model for subject: TCGA-BC-A10R ...
@ 2018-12-11 16:25:58: Building Bayesian Gaussian Mixture model for subject: TCGA-BC-A10S ...
@ 2018-12-11 16:26:07: Building Bayesian Gaussian Mixture model for subject: TCGA-BC-A10T ...
@ 2018-12-11 16:26:16: Building Bayesian Gaussian Mixture model for subject: TCGA-BC-A10U ...
@ 2018-12-11 16:26:24: Summerzie GMM models ...
@ 2018-12-11 16:26:24: Reports were saved into "summary_report.txt".
@ 2018-12-11 16:26:24: Writing to "TCGA-BC-A10Q.results.txt" ...
@ 2018-12-11 16:26:31: Writing to "TCGA-BC-A10R.results.txt" ...
@ 2018-12-11 16:26:38: Writing to "TCGA-BC-A10S.results.txt" ...
@ 2018-12-11 16:26:44: Writing to "TCGA-BC-A10T.results.txt" ...
@ 2018-12-11 16:26:50: Writing to "TCGA-BC-A10U.results.txt" ...

```                        
#### Output

```text
$ head TCGA-BC-A10Q.results.txt

#Prob_of_0: Probability of CpG belonging to un-methylation group
#Prob_of_1: Probability of CpG belonging to semi-methylation group
#Prob_of_2: Probability of CpG belonging to full-methylation group
#Assigned_lable: -1 = 'unsigned', 0 = 'un-methylation', 1 = 'semi-methylation', 2 = 'full-methylation'
Probe_ID	Beta_value	Prob_of_1	Prob_of_0	Prob_of_2	Assigned_lable
cg00000029	0.3469	1.0	5.8450778777948285e-31	1.894620863516962e-21	1
cg00000165	0.1517	0.997490700482822	0.0025092995171780517	4.9721760843021214e-39	-1
cg00000236	0.8479	0.08997964310682534	5.321555209091369e-221	0.9100203568931747	2
cg00000289	0.6658	0.9989955857078533	1.3294534692696416e-131	0.0010044142921467462	-1
cg00000292	0.6913	0.9932743887459475	1.2610051357010536e-142	0.006725611254052414	-1
```

####
Below histogram and piechart showed the proportion of CpGs assigned to "Un-methylated", 
"Semi-methylated" and "Full-methylated".

![trichotmize.png](https://github.com/liguowang/cpgtools/blob/master/img/trichotmize.png)

## <a name="p4"></a>Part 4: Comparison of differential CpG analysis tools

![Diff_analysis_table.png](https://github.com/liguowang/cpgtools/blob/master/img/Diff_analysis_table.png)
  

## <a name="p5"></a>Part 5: Contact Information
For questions, bugs, feature requests, and more. Created a ticket [here](https://github.com/liguowang/cpgtools/issues):
https://github.com/liguowang/cpgtools/issues
