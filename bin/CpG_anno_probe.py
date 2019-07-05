#!/usr/bin/env python
"""
Description
-----------
This program add comprehensive annotation information to each 450K/850K probe ID.
"""

import sys,os
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
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
	parser.add_option("-i","--input_file",action="store",type="string",dest="input_file",help="Input data file (Tab separated) with certain column containing 450K/850K array CpG IDs. This file can be regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-a","--annotation",action="store",type="string",dest="anno_file",help="Annotation file. This file can be regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
	parser.add_option("-p","--probe_column",action="store",type='int', dest="probe_col",default=0, help="The number specifying which column contains probe IDs. Note: the column index starts with 0. default=%default.")
	parser.add_option("-l", "--header", action="store_true", dest="header", default=False, help="Input data file has a header row.")
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
			print (l + '\t' +  '\t'.join(header), file=OUT)
		else:
			if options.probe_col >= len(f):
				print ("Error: column ID must be smaller than %d!" % len(f), file=sys.stderr)
				sys.exit(0)
			cgid = f[options.probe_col]
			try:
				print (l + '\t'  + data[cgid],file=OUT)
			except:
				print (l + '\t' + '\t'.join(['NA']*len(header)), file=OUT)
	OUT.close()		

if __name__=='__main__':
	main()
















		
