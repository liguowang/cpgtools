#!/usr/bin/env python
"""
Description
-----------
This program extracts annotation information for Illumina Epic array IDs. 
The first column of input file must contain CpG IDs (such as "cg18478105")
"""

import sys,os
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.1.8"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def read_annotation(infile):
	head = []
	cpg_infor = {}
	for l in ireader.reader(infile):
		if l.startswith('probeID'):
			head = l.split()[1:]
		else:
			f = l.split()
			cgid = f[0]
			anno = '\t'.join(f[1:])
			cpg_infor[cgid] = anno
	return (head, cpg_infor)
	
def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file with the first column containing CpG IDs. This file can be regular text file or compressed file (*.gz, *.bz2) or accessible url.")
	parser.add_option("-a","--annotation",action="store",type="string",dest="anno_file",help="Annotation file downloaded from \"https://github.com/liguowang/cpgtools/blob/master/epic/MethylationEPIC.tsv.gz\".")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
	parser.add_option("-l", "--header", action="store_true", dest="header", default=False, help="Data file has a header row.")
	(options,args)=parser.parse_args()
	
	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)

	if not (options.anno_file):
		print (__doc__)
		parser.print_help()
		sys.exit(102)
				
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)	

	if not os.path.isfile(options.input_file):
		print ("Input data file \"%s\" does not exist\n" % options.input_file) 
		sys.exit(104)
	if not os.path.isfile(options.anno_file):
		print ("Input annotation file \"%s\" does not exist\n" % options.input_file) 
		sys.exit(105)
			
	printlog("Read annotation file \"%s\" ..." % (options.anno_file))	
	(header, data)= read_annotation(options.anno_file)
	
	OUT = open(options.out_file + '.anno.txt','w')
	printlog("Add annotation information to \"%s\" ..." % (options.input_file))	
	line_num = 0
	for l in ireader.reader(options.input_file):
		line_num += 1
		f = l.split()
		if (line_num == 1 and options.header):
			first_col = f[0]
			other_cols = f[1:]
			print (first_col + '\t' + '\t'.join(header) + '\t' + '\t'.join(other_cols), file=OUT)
		else:
			cgid = f[0]
			other_cols = f[1:]
			try:
				print (cgid + '\t' + data[cgid] + '\t' + '\t'.join(other_cols),file=OUT)
			except:
				print (cgid + '\t' + '\t'.join(['NA']*len(header)) + '\t' + '\t'.join(other_cols),file=OUT)	
	OUT.close()		

if __name__=='__main__':
	main()
















		
