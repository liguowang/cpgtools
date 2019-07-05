#!/usr/bin/env python
"""
#=========================================================================================
This program picks the top N rows (according to standard deviation) from the input file.
The resulting file can be used for clustering/PCA analysis

Example of input data file
---------------------------
CpG_ID	Sample_01	Sample_02	Sample_03	Sample_04
cg_001	0.831035	0.878022	0.794427	0.880911
cg_002	0.249544	0.209949	0.234294	0.236680
cg_003	0.845065	0.843957	0.840184	0.824286
"""

import sys,os
import collections
import subprocess
import numpy as np
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED
import pandas as pd

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
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Tab separated data frame file containing beta values with the 1st row containing sample IDs and the 1st column containing CpG IDs.")
	parser.add_option("-c","--count",action="store",type='int', dest="cpg_count", default=1000, help="Number of most variable CpGs (ranked by standard deviation) to keep. default=%default" )
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
	(options,args)=parser.parse_args()
	
	print ()
	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)
	
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)	
	
	printlog("Reading input file: \"%s\"" % (options.input_file))
	df1 = pd.read_table(options.input_file, index_col = 0)
	
	#remove any rows with NAs
	df2 = df1.dropna(axis=0, how='any')
	printlog("%d rows with missing values were removed." % (len(df1) - len(df2)))
	
	#calculate stdev for each row
	row_stds = df2.std(axis=1)
	df2.loc[:, 'Stdev'] =  row_stds
	
	#sorted data frame by stdev (decreasingly). Then take the top count,. Then remove Stdev column
	printlog("Sorting by the standard deviation (decreasingly) ... ")
	df3 = df2.sort_values(by=['Stdev'], ascending=False)
	
	printlog("Data frame with sorted Stdev is saved to file: %s" % options.out_file + '.sortedStdev.tsv')
	df3.to_csv(options.out_file + '.sortedStdev.tsv', sep = "\t",float_format='%.6f')
		
	df4 = df3[0:options.cpg_count].drop('Stdev',axis=1)
	printlog("Top %d rows of Data frame is saved to file: %s" % (options.cpg_count, options.out_file + '.sortedStdev.topN.tsv'))
	df4.to_csv(options.out_file + '.sortedStdev.topN.tsv', sep="\t",float_format='%.6f')

	

if __name__=='__main__':
	main()	
