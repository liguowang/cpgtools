#!/usr/bin/env python
"""
#=========================================================================================
This program assign CpG to gene. Follows the "Basel plus extension" rule used by
GREAT (http://great.stanford.edu/public/html/index.php)
#=========================================================================================
"""


import sys,os
import collections
import subprocess
import numpy as np
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.1.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def main():
	print (__doc__)
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_files",help="Input CpG file(s) in BED format. BED file must have at least 3 columns: Chrom - the name of the chromosome or scaffold, chromStart - the starting position of CpG in the chromosome or scaffold, chromEnd - the ending position of CpG in the chromosome or scaffold. Note: the first base in a chromosome is numbered 0. BED file can be regular or compressed by 'gzip' or 'bz'. The barplot figures will NOT be generated if you provide more than 12 samples (bed files). [required]")
	parser.add_option("-r","--refgene",action="store",type="string",dest="gene_file",help="Reference gene model in standard BED-6 format (https://genome.ucsc.edu/FAQ/FAQformat.html#format1). ")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of output file. [required]")
	(options,args)=parser.parse_args()
	
	print ()

	if not (options.input_files):
		#print ('You must specify input file(s)',file=sys.stderr)
		parser.print_help()
		sys.exit(101)
	if not (options.gene_file):
		#print ('You must specify the chrom size file',file=sys.stderr)
		parser.print_help()
		sys.exit(102)
	if not (options.out_file):
		#print ('You must specify the output file',file=sys.stderr)
		parser.print_help()
		sys.exit(103)	
