import sys
from cpgmodule import ireader
import collections
from time import strftime
from bx.intervals import *
import numpy as np
from cpgmodule import ireader

def revcomp(dna):
	'''reverse complement DNA sequences'''
	tab = str.maketrans('ACGTNX*-','TGCANX*-')
	return dna.upper().translate(tab)[::-1]

def colors(n):
	'''
	return a list containing n colors
	'''
	if n >12 or n < 1:
		print("n must be in [1,12]", file=sys.stderr)
		return None
		
	color_12=['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']
	color_11=['#276419','#4d9221','#7fbc41','#b8e186','#e6f5d0','#f7f7f7','#fde0ef','#f1b6da','#de77ae','#c51b7d','#8e0152']
	color_10=['#276419','#4d9221','#7fbc41','#b8e186','#e6f5d0','#fde0ef','#f1b6da','#de77ae','#c51b7d','#8e0152']
	color_9 =['#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221']
	color_8 =['#c51b7d','#de77ae','#f1b6da','#fde0ef','#e6f5d0','#b8e186','#7fbc41','#4d9221']
	color_7 =['#c51b7d','#e9a3c9','#fde0ef','#f7f7f7','#e6f5d0','#a1d76a','#4d9221']
	color_6 =['#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221']
	color_5 =['#d01c8b','#f1b6da','#f7f7f7','#b8e186','#4dac26']
	color_4 =['#d01c8b','#f1b6da','#b8e186','#4dac26']
	color_3 =['#e9a3c9','#f7f7f7','#a1d76a']
	color_2 =['blue','red']
	color_1 =['blue']
	
	tmp=[color_1,color_2,color_3,color_4,color_5,color_6,color_7,color_8,color_9,color_10,color_11,color_12]
	return ["'" + i + "'" for i in tmp[n-1]]

def printlog (mesg):
	'''print progress message'''
	mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
	print (mesg, file=sys.stderr)
        
def chrom_count(infile):
	'''
	count chrom frequencies from BED file
	'''
	chrom_count = collections.defaultdict(int)
	for l in ireader.reader(infile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		if len(f)< 3:
			print ("BED has at lesat 3 columns. Skip: " + l, file=sys.stderr)
			continue
		try:
			start = int(f[1])
			end = int(f[2])
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue
		except:
			print ("Not in valid BED format. Skip:" + l, file=sys.stderr)
			continue
		
		chrom_count[f[0]] += 1
	return chrom_count

def read_chromSize(infile):
	'''
	read chromosome size file (tab/space separated plain text file).
	chr1    249250621
	chr2    243199373
	chr3    198022430
	chr4    191154276
	'''
	names = []
	sizes = []
	for l in ireader.reader(infile):
		if l.startswith('#'):
			continue
		f = l.split()
		if len(f) !=2:
			continue
		names.append(f[0])
		sizes.append(int(f[1]))
	return (names, sizes)

def equal_split(st, end, n):
	'''
	Equally split range(st,end) into n parts
	'''
	lst = []
	if end - st < n:
		return []
	stepSize = round((end - st)*1.0/n)
	count = 1
	
	a = st
	while count <= n:
		b = a + stepSize
		lst.append((a,b))
		a = b
		count += 1
	return lst


def read_CpG_bed(cpgfile):
	'''
	cpgfile: CpG BED file should have at least 3 columns (Chrom, chromStart, chromEnd).
	Note: chromEnd correspond to the genomic position methylated C.
	beta value is placed at the 4th column, if there is no 4th column (or the 4th column
	is not a number), beta set to 1.
	Additional columns are ignored. 
	'''
	cpg_ranges = {}
	for l in ireader.reader(cpgfile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		if len(f) < 3:
			print ("BED has at lesat 6 columns. Skip: " + l, file=sys.stderr)
			continue
		
		chrom = f[0]
		start = int(f[1])
		end = int(f[2])
		if start > end:
			print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
			continue		

		try:
			beta = float(f[4])
		except:
			beta = 1.0
		try:
			strand = f[5]
		except:
			strand = '+'		

		if chrom not in cpg_ranges:
			cpg_ranges[chrom] = IntervalTree()
		if strand == '+':
			cpg_ranges[chrom].insert_interval( Interval(start, end, value=beta))
		elif  strand == '-':
			cpg_ranges[chrom].insert_interval( Interval(end, end+1, value=beta))
		
	return cpg_ranges


def read_region_bed(bedfile):
	'''
	bedfile file should have at least 3 columns (Chrom, chromStart, chromEnd).
	if no strand information found in the 6th column. All regions will be 
	considered on "+" strand. 
	'''
	for l in ireader.reader(bedfile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		
		try:
			chrom = f[0]
			start = int(f[1])
			end = int(f[2])
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue		
		except:
			print ("BED has at lesat 3 columns. Skip: " + l, file=sys.stderr)
		try:
			strand = f[5]
		except:
			strand = "+"
		
		yield(chrom, start, end, strand)

def read_bed_as_list(bedfile):
	'''
	bedfile file should have at least 3 columns (Chrom, chromStart, chromEnd).
	if no strand information found in the 6th column. All regions will be 
	considered on "+" strand. 
	'''
	lst = []
	for l in ireader.reader(bedfile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		
		try:
			chrom = f[0]
			start = int(f[1])
			end = int(f[2])
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue		
		except:
			print ("BED has at lesat 3 columns. Skip: " + l, file=sys.stderr)		
		lst.append([chrom, start, end])
	return lst

def coverage_over_range(lst, cpg_ranges):
	'''
	Calculate relative methylation density
	lst = list of (chr,start,end, strand)
	cpg_ranges is returned by read_CpG_bed
	'''
	
	results = collections.defaultdict(list)
	beta_signals = {}
	for (chr,st,end, strand) in lst:
		if chr not in cpg_ranges:
			continue
		
		span = end - st
		tmp = cpg_ranges[chr].find(st, end)	#eg: [Interval(3, 40, value=3), Interval(13, 50, value=4)]
		for i in tmp:
			if strand == '+':
				CpG_to_origin = round((i.end - (st+1))*100/span)
			if strand == '-':
				CpG_to_origin = abs(round((i.end - end)*100/span))
			CpG_beta = i.value
			results[CpG_to_origin].append(CpG_beta)
	for k,v in results.items():
		beta_signals[k] = round(np.mean(v),4)
	return beta_signals


def count_over_range(lst, cpg_ranges):
	'''
	Calculate how many CpGs are located in lst
	lst = list of (chr,start,end)
	cpg_ranges is returned by read_CpG_bed
	'''
	
	total_size = 0	#total nucleotides of list of genomic regions
	total_count = 0 #total CpGs in list of genomic regions
	for (chr,st,end) in lst:
		total_size += (end - st)
		if chr not in cpg_ranges:
			continue
		tmp = cpg_ranges[chr].find(st, end)	#eg: [Interval(3, 40, value=3), Interval(13, 50, value=4)]
		total_count += len(tmp)
	return(total_size,total_count)

def read_grp_file1(gfile):
	'''
	read group file. Group file define the biological groups of data matrix file. 
	(1) It must has header
	(2) It must have two columns:
		* 1st column: sample names. samples names should be unique, and they must be exactly the same as the first row of beta matrix file.
		* 2nd column: group IDs. 
	(3) columns must be separated by ","
	
	For example:
	
	sampleID,groupID
	Normal_1,1
	Normal_2,1
	Normal_3,1
	Tumor_1,2
	Tumor_2,2
	Tumor_3,2
	'''
	samples = []
	groups = []
	line_num = 0
	for l in ireader.reader(gfile):
		l = l.replace(' ','')
		line_num += 1
		f = l.split(',')
		if len(f) < 2:
			print ("Group fle must have 2 columns!", file=sys.stderr)
			sys.exit(1)
		if line_num == 1:
			continue
		else:		
			samples.append(f[0])
			groups.append(f[1])
	
	tmp = collections.Counter(samples)
	if tmp.most_common(1)[0][1] > 1:
		print ("Sample names are not unique!", file=sys.stderr)
		sys.exit(0)
		
	return(samples, groups)

def read_grp_file2(gfile):
	'''
	read group file. Group file define the biological groups of data matrix file. 
	(1) It must has header
	(2) It must have at least two columns:
		* 1st column: sample names. samples names should be unique, and they must be exactly the same as the first row of beta matrix file.
		* 2nd column: group IDs. 
		* additional columns can be included to indicate co-variables. 
	(3) columns must be separated by ","
	
	For example:
	
	sampleID,survival,Sex
	Normal_1,1,1
	Normal_2,1,2
	Normal_3,1,1	
	Tumor_1,2,1
	Tumor_2,2,2
	Tumor_3,2,1
	...
	...
	'''
	samples = []
	covar_values = []
	covar_names = []
	covars = collections.defaultdict(dict)
	line_num = 0
	
	covar_values = collections.defaultdict(list)	#continue variable or categorical variable. key is name, valu list of values
	cutoff = 0.5	#ratio of number of unique values to the total number of unique values
	for l in ireader.reader(gfile):
		l = l.replace(' ','')
		line_num += 1
		f = l.split(',')
		if len(f) < 2:
			print ("Group fle has at lesat 2 columns!", file=sys.stderr)
			sys.exit(1)
		if line_num == 1:
			covar_names = f[1:]
		else:
			sample_id = f[0]
			samples.append(sample_id)
			row_values = f[1:]
			
			for a,b in zip(covar_names, row_values):
				covars[a][sample_id] = b
				covar_values[a].append(b)
		
	tmp = collections.Counter(samples)
	if tmp.most_common(1)[0][1] > 1:
		print ("Sample names are not unique!", file=sys.stderr)
		sys.exit(0)
	
	#tell if a covariable is continuous or categorical
	covar_types = {}
	for k,v in covar_values.items():
		if ( 1.0*len(set(v)) / len(v) ) > cutoff:
			covar_types[k] = 'continuous'
		else:
			covar_types[k] = 'categorical'
		
	return(samples, covar_names, covars, covar_types)
			
def stats_over_range(cpg_ranges, chrom, st, end):
	'''
	Basic statistics about range
	'''
	
	stats = []

	if chrom not in cpg_ranges:
		return ['NA']*6
	
	tmp = []
	overlaps = cpg_ranges[chrom].find(st, end)
	for i in overlaps:
		tmp.append(i.value)
			
	if len(tmp) == 0:
		return ['NA']*6

	try:
		i_count = len(overlaps)
	except:
		i_count = 'NA'
	
	try:
		i_min = round(min(tmp),4)
	except:
		i_min = 'NA'
	
	try:
		i_max = round(max(tmp),4)
	except:
		i_max = 'NA'
	
	try:
		i_mean = round(np.mean(tmp),4)
	except:
		i_mean = 'NA'
	
	try:
		i_median = round(np.median(tmp),4)
	except:
		i_median = 'NA'
	
	try:
		if len(tmp) > 1:
			i_std = round(np.std(tmp, ddof=1),4)
		else:
			i_std = 'NA'
	except:
		i_std = 'NA'
	
	return [i_count, i_min, i_max, i_mean, i_median, i_std]

def load_pickle_obj():
    with open('./id2chr.pkl', 'rb') as f:
        return pickle.load(f)


"""			
def read_CpG_bed(cpgfile,genefile, bin_count = 100):
	'''
	cpgfile: CpG BED file (at least 3 columns).
	genefile: gene BED file (at least 6 columns, must have strand information). 
	'''
	cpg_ranges = {}
	for l in ireader.reader(cpgfile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		if len(f)< 3:
			print ("BED has at lesat 3 columns. Skip: " + l, file=sys.stderr)
			continue
		try:
			chrom = f[0]
			start = int(f[1])
			end = int(f[2])
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue
		except:
			print ("Not in valid BED format. Skip:" + l, file=sys.stderr)
			continue

		if chrom not in cpg_anges:
			cpg_ranges[chrom] = IntervalTree()
		cpg_ranges[chrom].insert_interval( Interval( int(start), int(end)))
		
	#return cpg_ranges
	
	cpg_profile = []	#list of list = [CpG_count across bins]
	for l in ireader.reader(genefile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		if len(f)< 6:
			print ("Gene BED has at lesat 6 columns. Skip: " + l, file=sys.stderr)
			continue
		try:
			chrom = f[0]
			tss_start = int(f[1])
			tss_end = int(f[2])
			strand = f[5]
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue
		except:
			print ("Not in valid BED format. Skip:" + l, file=sys.stderr)
			continue
		
		#
		if chrom not in cpg_ranges:
			continue
			
		genomic_size = tss_end - tss_start
		window_start = tss_start - int(genomic_size/2.0)	#extend upstream half gene size
		window_end = tss_end + int(genomic_size/2.0)		#extend downstream half gene size
		if window_start < 0:
			window_start = 0
		
		
		bins = equal_split(widow_start, window_end, bin_count)
		if len(bins) == 0: continue
		
		cpg_counts = []	#CcG count in each bin
		for (bin_st, bin_end) in bins:
			tmp = cpg_ranges[chrom].find(bin_st, bin_end)
			cpg_counts.append(len(tmp))
		
		if strand == '-':
			cpg_counts = cpg_counts[::-1]
		cpg_profile.append(cpg_counts)
	
	return np.array(cpg_profile).means(axis=0)	
"""		
	
