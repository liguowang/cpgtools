#!/usr/bin/env python
"""
Description
-----------
This program calculates the distribution of CpG over user-specified genomic regions. 

Notes
------
 1. A maximum of 10 BED files (define 10 different genomic regions) can be analyzed
    together. 
 2. The *order* of BED files is important (i.e. considered as "priority order"). Overlapped
    genomic regions will be kept in the BED file with the highest priority and removed
    from BED files of lower priorities.  For example, users provided 3 BED files via  "-i
    promoters.bed,enhancers.bed,intergenic.bed", then if an enhancer region is overlapped
    with promoters, *the overlapped part* will be removed from "enhancers.bed".
 3. BED files can be regular or compressed by 'gzip' or 'bz'.
#=========================================================================================
"""


import sys,os
import collections
import subprocess
import numpy as np
#import re
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
	parser.add_option("-i","--cpg",action="store",type="string",dest="cpg_file",help="BED file specifying the C position. This BED file should have at least 3 columns (Chrom, ChromStart, ChromeEnd).  Note: the first base in a chromosome is numbered 0. This file can be a regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-b","--bed",action="store",type="string",dest="bed_files",help="List of BED files specifying the genomic regions.")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
	(options,args)=parser.parse_args()
	
	print ()

	if not (options.cpg_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)
		
	if not (options.bed_files):
		print (__doc__)
		parser.print_help()
		sys.exit(101)
				
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(102)	

	FOUT = open(options.out_file + '.txt','w')
	ROUT = open(options.out_file + '.r','w')
	
	#step1: read CpG file
	printlog("Reading CpG file: \"%s\"" % (options.cpg_file))
	cpg_ranges = read_CpG_bed(options.cpg_file)
	
	#step2: check BED file
	printlog("Checking BED files: \"%s\"" % (options.bed_files))
	input_bed_files = options.bed_files.replace(' ','').split(',')
	for i in input_bed_files:
		if os.path.exists(i):
			print("\t%s" % i, file=sys.stderr)
		else:
			print("\"%s\" does not exist!" % i, file=sys.stderr)
			sys.exit(103)
	
	#step3: read, merge, and subtract BED file
	dat = {}
	result = [("Priority_order", "Name", "Number_of_regions", "Size_of_regions(bp)", "CpG_raw_count", "CpG_count_per_KB")]
	
	#step3.1: read the first BED file
	i = 0
	printlog("Reading BED file: \"%s\"" % (input_bed_files[i]))
	file_name = os.path.basename(input_bed_files[i])
	tmp = read_bed_as_list(input_bed_files[i])
	printlog("Merging overlap entries in BED file: \"%s\"" % (input_bed_files[i]))
	dat[i] = BED.unionBed3(tmp)
	printlog("Counting CpGs ...")
	(size,count) = count_over_range(dat[i], cpg_ranges)
	result.append([str(i), file_name, len(dat[i]), size, count, count*1000.0/size])	#Class, number_of_region, size_of_region, CpG_raw_count, CpG_count_perKb
	
	#step3.2: read the remaining BED files
	for i in range(1, len(input_bed_files)):
		printlog("Reading BED file: \"%s\"" % (input_bed_files[i]))
		file_name = os.path.basename(input_bed_files[i])
		tmp = read_bed_as_list(input_bed_files[i])
		printlog("Merging overlap entries in BED file: \"%s\"" % (input_bed_files[i]))
		dat[i] = BED.unionBed3(tmp)
		
		for j in range(0,i):
			printlog("Subtract \"%s\" from \"%s\"" % (input_bed_files[j], input_bed_files[i]))
			dat[i] = BED.subtractBed3(dat[i],  dat[j])
		(size,count) = count_over_range(dat[i], cpg_ranges)
		result.append([str(i), file_name, len(dat[i]), size, count, count*1000.0/size])
			
	print('\n')
	names=[]	#[0,1,2,3,4,...]
	labels = []	#[bed names]
	density=[]
	for tmp in result:
		print ('\t'.join([str(i) for i in tmp]), file=FOUT)
		names.append(tmp[0])
		labels.append(tmp[1])
		density.append(tmp[5])
	FOUT.close()
	
	print("name = c(%s)" % ','.join(['"' + i + '"' for i in names[1:]]), file=ROUT)
	print("values = c(%s)" % ','.join([str(i) for i in density[1:]]), file=ROUT)
	print ('pdf("%s", width=8, height=6)' % (options.out_file + '.pdf'), file=ROUT)
	print ('layout(matrix(c(1,1,2,1,1,2), nrow=2, byrow=TRUE))', file=ROUT)
	print ('barplot(values,names.arg=name,col=c(%s),ylab="CpG per Kb")' % ','.join(colors(len(input_bed_files))), file=ROUT)
	print ("plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')", file=ROUT)
	for name,label in zip(names[1:], labels[1:]):
		x_pos = 0.0
		y_pos = 1-(int(name)*20.0 +5)/200 
		print ("text(x=%f, y=%f, labels=c(\"%s = %s\"),adj=c(0,0))" % (x_pos, y_pos,name,label), file=ROUT)
	print ('dev.off()', file=ROUT)
	
	ROUT.close()
	
	printlog("Running R script ...")
	try:
		subprocess.call("Rscript " + options.out_file + '.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
		pass		

if __name__=='__main__':
	main()	
	
