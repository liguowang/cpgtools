"""
Aggregate proportion values of a list of CpGs that located in give genomic regions
(eg. CpG islands, promoters, exons, etc).

Outlier CpG will be removed if the probability of observing its proportion vlaue is less
than p-cutoff. For example, if alpha set to 0.05 and there are 10 CpGs (n = 10) located in a 
particular genomic region, the p-cutoff of this genomic region is 0.005 (0.05/10). Supposing
the total reads mapped to this region is 100, out of which 25 are methylated reads (i.e
regional methylation level (beta) = 25/100 = 0.25) 

The probability of observing CpG (3,10) is :
	pbinom(q=3, size=10, prob=0.25) = 0.7759
The probability of observing CpG (0,10) is :
	pbinom(q=0, size=10, prob=0.25) = 0.05631
The probability of observing CpG (16,21) is :
	pbinom(q=16, size=21, prob=0.25, lower.tail=FALSE) = 1.19e-07 (outlier)


**Example of input file**

Chrom	Start	End	score
chr10	100017748	100017749	3,10	
chr10	100017769	100017770	0,10	
chr10	100017853	100017854	16,21

"""

import sys,os
import collections
import subprocess
import numpy as np
from scipy.stats import binom

from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED
import pandas as pd
from bx.intervals import *

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"



def buildIntervalTree(bed_file):
	'''
	Build interval tree from BED file. Input BED file must have at least 4 columns
	'''        
	ranges={}
	printlog("reading "+ bed_file + '...')
	for line in ireader.reader(bed_file):
		if line.startswith("track"):continue
		if line.startswith("#"):continue
		if line.startswith('browser'):continue   
		if line.startswith('Chrom'):continue	
		fields = line.rstrip('\n ').split()
		if len(fields) < 4:
			continue
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])
		score = fields[3]
		
		if start < 0:
			continue
		if end < 0:
			continue
		if start > end:
			continue
		

		if chrom not in ranges:
			ranges[chrom] = Intersecter()
			ranges[chrom].add_interval( Interval( start, end, value=score) )
		else:
			ranges[chrom].add_interval( Interval( start, end, value=score) )
	return ranges

def findIntervals(chrom, start, end, obj, a = 0.01):
	'''
	obj is the IntervalTree object returned by "buildIntervalTree.
	'''
	hits = []	# list of proportion values
	
	if chrom not in obj:
		return hits
	else:
		overlaps = obj[chrom].find(int(start), int(end))
		for i in overlaps:
			hits.append(i.value)
	if len(hits) == 0:
		return(['N/A']*6)
	
	
	methyl = []	#list of methylated read for each CpG
	total = []	#list of total read for each CpG

	for h in hits:
		m, t = h.split(',')
		methyl.append(int(m))
		total.append(int(t))
	ori_CpG_count = len(total)		#number of CpGs of a region
	p_cut = a / ori_CpG_count
	ori_methyl_sum = np.sum(methyl)	#total reads of a region
	ori_total_sum = np.sum(total)	#total methylated reads of a region
	
	if ori_total_sum == 0:
		return(['N/A']*6)
	if ori_methyl_sum == 0 or ori_methyl_sum == ori_total_sum:
		return([ori_CpG_count, ori_methyl_sum, ori_total_sum, ori_CpG_count, ori_methyl_sum, ori_total_sum])
	
	
	region_beta = ori_methyl_sum/ori_total_sum	#average methylation level of *region*, equivalent to prob in binomial 

		
	new_methyl = []
	new_total = []
	for m, t in zip(methyl, total):
		p = binom.cdf(k = m, n = t, p = region_beta)
		#print (p, m, t)
		if p < p_cut:
			continue
		if (1.0 - p) < p_cut:
			continue
		new_methyl.append(m)
		new_total.append(t)
	new_CpG_count = len(new_total)
	new_methyl_sum = np.sum(new_methyl)
	new_total_sum = np.sum(new_total)
	
	return([new_CpG_count, new_methyl_sum, new_total_sum, ori_CpG_count, ori_methyl_sum, ori_total_sum])
	
	

def main():
	
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input",action="store",type="string",dest="input_file",help="Input CpG file in BED format. The first 3 columns contain \"Chrom\", \"Start\", and \"End\". The 4th column contains proportion values.")
	parser.add_option("-a","--alpha",action="store",type='float', dest="alpha_cut", default=0.05, help="The chance of mistakingly assign a particular CpG as an outlier for each genomic region. default=%default" )
	parser.add_option("-b","--bed",action="store",type="string",dest="bed_file",help="BED3+ file specifying the genomic regions.")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
	(options,args)=parser.parse_args()

	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)
	
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)
	
	if options.alpha_cut < 0:
		options.alpha_cut = 0.05
	if options.alpha_cut > 1:
		options.alpha_cut = 1	

	tree = buildIntervalTree(options.input_file)
	
	OUT = open(options.out_file,'w')
	print ("#chrom\tstart\tend\tN_CpG_filtered\tN_methyl_filtered\tN_total_filtered\tN_CpG_ori\tN_methy_ori\tN_total_ori", file=OUT)	
	for line in ireader.reader(options.bed_file):
		line = line.strip()
		if line.startswith("track"):continue
		if line.startswith("#"):continue
		if line.startswith('browser'):continue   
		if line.startswith('Chrom'):continue	

		
		f = line.split()
		if len(f) < 3:
			continue
		try:
			chrom = f[0]
			start = int(f[1])
			end = int(f[2])
		except:
			continue
		
		a = findIntervals(chrom, start, end, tree, a = options.alpha_cut)
		if len(a) == 0:
			print ('\t'.join(f[0:3]) + '\t' + '\t'.join( ['N/A']*6), file=OUT)
		else:
			print ('\t'.join(f[0:3]) + '\t' + '\t'.join([str(i) for i in a]), file=OUT)
	OUT.close()
	
if __name__=='__main__':
	main()