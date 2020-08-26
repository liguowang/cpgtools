#!python

"""
Description
-----------
This program annotates CpG by its position.

Notes:
- Input CpG and BED files must have at least three columns
- If multiple regions from the annotation BED file are overlapped with the **same**
  CpG site, their names will be concatenated together. 
   
"""

import sys,os
import collections
import subprocess
import numpy as np
from os.path import basename
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
__version__="0.1.8"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"



def buildIntervalTree(bed_file, window_size = 0):
	'''
	Build interval tree from annotation BED file.
	window : add this to the middle of each region.
	'''        
	ranges={}
	printlog("Build interval tree from annotation file: %s ..." % bed_file)
	for line in ireader.reader(bed_file):
		if line.startswith("track"):continue
		if line.startswith("#"):continue
		if line.startswith('browser'):continue   
		fields = line.rstrip('\n ').split()
		if len(fields) < 3:
			continue
		chrom = fields[0]
		start = int(fields[1])
		end = int(fields[2])
		#middle position
		mid = int(start + (end - start)/2.0 )
		if start < 0:
			continue
		if end < 0:
			continue
		if start > end:
			continue
		
		# window start position
		extension = int(window_size * 0.5)
		w_start = mid - extension
		if w_start < start:
			w_start = start
		
		# window end position
		w_end = mid + extension
		if w_end > end:
			w_end = end
			
		if len(fields) >= 4:
			name = fields[3]
		else:
			name = fields[0] + ':' + fields[1] + '-' + fields[2]

		if chrom not in ranges:
			ranges[chrom] = Intersecter()
			ranges[chrom].add_interval( Interval( start, end, value=name) )
		else:
			ranges[chrom].add_interval( Interval( start, end, value=name) )
	return ranges

def findIntervals(chrom, start, end, obj):
	'''
	obj is the IntervalTree object returned by "buildIntervalTree.
	'''
	hits = set()
	if chrom not in obj:
		return hits
	else:
		overlaps = obj[chrom].find(int(start), int(end))
		for i in overlaps:
			hits.add(i.value)
	return sorted(hits)

def main():
	
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input_file",action="store",type="string",dest="input_file",help="Input CpG file in BED3+ format.")
	parser.add_option("-a","--annotation",action="store",type="string",dest="anno_file",help="Input annotation file in BED3+ format.")
	parser.add_option("-w","--window",action="store",type='int', dest="window_size", default=100, help="Size of window centering on the middle-point of each genomic region defined in the annotation BED file (i.e., window_size*0.5 will be extended to up- and down-stream from the middle point of each genomic region). default=%default" )
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="The prefix of the output file.")
	parser.add_option("-l", "--header", action="store_true", dest="header", default=False, help="If True, the first row of input CpG file is header. default=%default")
	(options,args)=parser.parse_args()

	
	if not (options.input_file):
		print (__doc__)
		#print ('You must specify input file(s)',file=sys.stderr)
		parser.print_help()
		sys.exit(101)
	if not (options.out_file):
		print (__doc__)
		#print ('You must specify the output file',file=sys.stderr)
		parser.print_help()
		sys.exit(102)	
	if not (options.anno_file):
		print (__doc__)
		#print ('You must specify the annotation file',file=sys.stderr)
		parser.print_help()
		sys.exit(103)
	tree = buildIntervalTree(options.anno_file, window_size = options.window_size)
	
	OUT = open(options.out_file + '.anno.txt','w')
	line_num = 0
	printlog("Reading CpG file: %s ..." % options.input_file)
	for line in ireader.reader(options.input_file):	 
		fields = line.rstrip('\n ').split()
		if len(fields) < 3:
			continue
		line_num += 1
		f = line.split()
		if (line_num == 1 and options.header):
			print (line + '\t' +  basename(options.anno_file), file=OUT)
		else:
			chrom = f[0]
			start = int(f[1])
			end = int(f[2])
			overlaps = findIntervals(chrom, start, end, tree)
			if len(overlaps) > 0:
				print (line + '\t' + ','.join(overlaps), file=OUT)
			else:
				print (line + '\tN/A', file=OUT)
	
	OUT.close()
		
if __name__=='__main__':
	main()
	