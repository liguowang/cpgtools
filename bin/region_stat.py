#!/usr/bin/env python
"""
#=========================================================================================
This program gives basic statistics for each genomic region. Add 6 columns to the BED file:
1. Number of CpGs detected in the genomic region
2. Min methylation level
3. Max methylation level
4. Average methylation level across all CpGs
5. Median methylation level across all CpGs
6. Standard deviation
#=========================================================================================
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
__version__="0.1.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="BED6 file specifying the C position. This BED file should have at least 6 columns (Chrom, ChromStart, ChromeEnd, Name, Beta_value, Strand).  Note: the first base in a chromosome is numbered 0. This file can be a regular text file or compressed file (*.gz, *.bz2) or accessible url.")
	parser.add_option("-r","--region",action="store",type="string",dest="region_file",help="BED file specificy genomic regions. This BED file should have at least 3 columns (Chrom, ChromStart, ChromeEnd).")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of output file.")
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
				
