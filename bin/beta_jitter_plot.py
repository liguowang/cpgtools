#!/usr/bin/env python
"""
Description
-----------
This program generates jitter plot (a.k.a. strip chart) and bean plot for each sample (column). 

Notes
-----
User must install the "beanplot" R library:
https://cran.r-project.org/web/packages/beanplot/index.html

Example of input
-----------------
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
	parser.add_option("-i","--input",action="store",type="string",dest="input_file",help="Tab separated data frame file containing beta values with the 1st row containing sample IDs and the 1st column containing CpG IDs.")
	parser.add_option("-f","--fraction",action="store",type='float', dest="fraction", default=0.5, help="Fraction of total data points (CpGs) used to generate jitter plot. Decrease this number if the jitter plot is over-crowded. default=%default" )
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
	if options.fraction < 0:
		options.fraction = 0.0
	if options.fraction > 1:
		options.fraction = 1.0	
	
	ROUT = open(options.out_file + '.r','w')
	print ('library(beanplot)', file = ROUT)
	print ('pdf(file=\"%s\", width=10, height=8)' % (options.out_file + '.pdf'),file=ROUT)
	
	
	
	printlog("Reading beta file: \"%s\"" % (options.input_file))
	df = pd.read_table(options.input_file)	
	samples = df.columns[1:]
	if options.fraction < 1.0:
		subset_file = options.out_file + '.sample.tsv'
		printlog("Sampling subset (%.2f%%) from file: \"%s\"" % (options.fraction*100, options.input_file))
		subset = df.sample(frac=options.fraction, replace=False, random_state=999)
		printlog("Saving subset (%.2f%%) to file: \"%s\"" % (options.fraction*100, options.out_file + '.sample.tsv'))
		subset.to_csv(options.out_file + '.sample.tsv', sep = "\t", index = False)
		
		print ('d = read.table(file="%s", sep="\\t", header=TRUE)' % (options.out_file + '.sample.tsv'), file=ROUT)
		print ('ll = list(%s)' % (','.join(['"' + s + '"' + ' = d$' + s for s in samples])), file=ROUT)
		print ('stripchart(ll,cex=0.1,col="#abd9e9", vertical=T, method=c("jitter"), ylab="Beta value",las=2, jitter=0.3,cex.names = 0.8)', file=ROUT)
		print ('beanplot(ll, cutmin=0,cutmax=1, border="#d01c8b",what=c(1,1,1,0),col=c(), las=2, add = TRUE)', file=ROUT)
	else:
		printlog("Using all data points in file: \"%s\"" % (options.input_file))
		print ('d = read.table(file="%s", sep="\\t", header=TRUE)' % (options.input_file), file=ROUT)
		print ('ll = list(%s)' % (','.join(['"' + s + '"' + ' = d$' + s for s in samples])), file=ROUT)
		print ('stripchart(ll,cex=0.1,col="#abd9e9", vertical=T, method=c("jitter"), ylab="Beta value", las=2, jitter=0.3,cex.names = 0.8)', file=ROUT)
		print ('beanplot(ll, cutmin=0,cutmax=1, border="#d01c8b",what=c(1,1,1,0),col=c(), las=2, add = TRUE)', file=ROUT)
	ROUT.close()
	
	try:
		subprocess.call("Rscript " + options.out_file + '.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
		pass        	
	
if __name__=='__main__':
	main()	
	
