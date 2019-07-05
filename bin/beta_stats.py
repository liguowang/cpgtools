#!/usr/bin/env python
"""
Description
-----------
This program gives basic information of CpGs located in each genomic region. 
It adds 6 columns to the input BED file:
 1. Number of CpGs detected in the genomic region
 2. Min methylation level
 3. Max methylation level
 4. Average methylation level across all CpGs
 5. Median methylation level across all CpGs
 6. Standard deviation
 
Example of input BED6+ file
---------------------------
chr22   44021512        44021513        cg24055475      0.9231  -
chr13   111568382       111568383       cg06540715      0.1071  +
chr20   44033594        44033595        cg21482942      0.6122  -
"""


import sys,os
import collections
import subprocess
import numpy as np
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="BED6+ file specifying the C position. This BED file should have at least 6 columns (Chrom, ChromStart, ChromeEnd, Name, Beta_value, Strand).  Note: the first base in a chromosome is numbered 0. This file can be a regular text file or compressed file (.gz, .bz2)")
	parser.add_option("-r","--region",action="store",type="string",dest="region_file",help="BED3+ file of genomic regions. This BED file should have at least 3 columns (Chrom, ChromStart, ChromeEnd).")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
	(options,args)=parser.parse_args()
	
	print ()

	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)

	if not (options.region_file):
		print (__doc__)
		parser.print_help()
		sys.exit(102)
				
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)	
	
	FOUT = open(options.out_file + '.txt','w')
	
	#step1: read CpG file
	printlog("Reading CpG file: \"%s\"" % (options.input_file))
	cpg_ranges = read_CpG_bed(options.input_file)
		
	#step2: read region file
	printlog("Reading BED file: \"%s\"" % (options.region_file))
	
	printlog("Writing to: \"%s\"" % (options.out_file + '.txt'))
	region_list = []
	for l in ireader.reader(options.region_file):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		if len(f) < 3:
			continue
		try:
			chrom = f[0]
			st = int(f[1])
			end = int(f[2])
		except:
			print (l + '\t' + '\t'.join(['NA']*6, file=FOUT))
			continue
		tmp = stats_over_range(cpg_ranges, chrom, st, end)
		print (l + '\t' + '\t'.join([str(i) for i in tmp]), file=FOUT)		
	
	FOUT.close()

if __name__=='__main__':
	main()	
				
