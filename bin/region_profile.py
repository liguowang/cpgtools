#!/usr/bin/env python
"""
#=========================================================================================
This program calculates the overall methylation level (i.e. average beta value) across 
particular genomic regions (eg. promoters, TF bindings). 
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
	parser.add_option("-r","--region",action="store",type="string",dest="region_file",help="BED file specificy genomic regions. This BED file should have at least 3 columns (Chrom, ChromStart, ChromeEnd). If the 6-th column does not exist, all regions will be considered as on \"+\" strand. ")
	parser.add_option("-d","--downstream",action="store",type="int",dest="downstream_size",default=2000,help="Size of extension to downstream. default=%default (bp)")
	parser.add_option("-u","--upstream",action="store",type="int",dest="upstream_size",default=2000,help="Size of extension to upstream. default=%default (bp)")
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
	ROUT = open(options.out_file + '.r','w')
	print ("\t".join(["Group","Relative_position(5'->3')", "Average_beta"]), file=FOUT)
	
	#step1: read CpG file
	printlog("Reading CpG file: \"%s\"" % (options.input_file))
	cpg_ranges = read_CpG_bed(options.input_file)
		
	#step2: read region file
	printlog("Reading BED file: \"%s\"" % (options.region_file))
	
	region_list = []
	for chrom, st, end, strand in read_region_bed(options.region_file):
		region_list.append((chrom, st, end, strand))
	region_list = list(set(region_list))

	printlog("Calculate average beta ...")
	s = coverage_over_range(region_list,cpg_ranges)
	for i in sorted(s):
		print ('\t'.join(["User_region", str(i), str(s[i])]), file=FOUT)	
	print ('User_region <- c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)
	user_region_datapoints = len(s)
	
	if options.upstream_size > 0:
		printlog("Get upstream regions of  \"%s\"" % (options.region_file))
		upstream_region = []
		for (chrom, st, end, strand) in region_list:
			if strand == '+':
				upstream_st = max(st - options.upstream_size, 0)
				upstream_end = st
				upstream_region.append((chrom, upstream_st, upstream_end, strand))
			elif strand == '-':
				upstream_st = end
				upstream_end = end + options.upstream_size
				upstream_region.append((chrom, upstream_st, upstream_end, strand))
		upstream_region = list(set(upstream_region))

	s = coverage_over_range(upstream_region,cpg_ranges)
	for i in sorted(s):
		print ('\t'.join(["Upstream_region", str(i), str(s[i])]), file=FOUT)	
	print ('Upstream_region <- c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)
	upstream_datapoints = len(s)


	if options.downstream_size > 0:
		printlog("Get downstream regions of  \"%s\"" % (options.region_file))
		downstream_region = []
		for (chrom, st, end, strand) in region_list:
			if strand == '+':
				downstream_st = end
				downstream_end = end + options.downstream_size
				downstream_region.append((chrom, downstream_st, downstream_end, strand))
			elif strand == '-':
				downstream_st = st
				downstream_end = max(st - options.downstream_size, 0)
				downstream_region.append((chrom, downstream_st, downstream_end, strand))
		downstream_region = list(set(downstream_region))
	s = coverage_over_range(downstream_region,cpg_ranges)
	for i in sorted(s):
		print ('\t'.join(["Downstream_region", str(i), str(s[i])]), file=FOUT)	
	print ('Downstream_region <- c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)
	downstream_datapoints = len(s)
	
	total_datapoints = upstream_datapoints + downstream_datapoints + user_region_datapoints
	print('\n')
	print ('pdf(file=\"%s\", width=6, height=6)' % (options.out_file + '.pdf'),file=ROUT)
	print ('plot(0:%d, c(Upstream_region, User_region, Downstream_region),ylim=c(0,1), xaxt="n",xlab="", ylab="Average methylation", type="l", col="red")' % (total_datapoints -1), file=ROUT)
	print ('abline(v = c(%d,%d),col="blue", lty="dashed")' % (upstream_datapoints-1, upstream_datapoints + user_region_datapoints - 1), file=ROUT)
	print ('abline(h = 0.5,col="grey", lty="dashed")', file=ROUT)
	print ('text(x=c(%d, %d), y=0.9, cex=0.7, labels=c("Upstream\\n(5\'->3\')","Downstream\\n(5\'->3\')"))' % (50, total_datapoints - 50), file=ROUT)
	print ('dev.off()',file=ROUT)	
	
	FOUT.close()
	ROUT.close()
	try:
		subprocess.call("Rscript " + options.out_file + '.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
		pass        
		
		

if __name__=='__main__':
	main()	
