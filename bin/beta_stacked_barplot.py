#!/usr/bin/env python
"""
Description
-----------
This program creates stacked barplot for each sample. The stacked barplot showing
the proportions of CpGs whose beta values are falling into these 4 ranges:
 * [0.00,  0.25] 	#first quantile
 * [0.25,  0.50]	#second quantile
 * [0.50,  0.75]	#third quantile
 * [0.75,  1.00]	#forth quantile

Example of input data file
---------------------------
CpG_ID	Sample_01	Sample_02	Sample_03	Sample_04
cg_001	0.831035	0.878022	0.794427	0.880911
cg_002	0.249544	0.209949	0.234294	0.236680
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
import pandas as pd

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def quarter_count(lst):
	"""
	count number of beta values falling into each quarter
	Note: beta >= 0 and beta <=1
	"""
	q1 = 0
	q2 = 0
	q3 = 0
	q4 = 0
	for i in lst:
		try:
			j = float(i)
		except:
			continue
		if not isinstance(j, float):
			continue
		if j < 0:
			continue
		elif j <= 0.25:
			q1 += 1
		elif j <= 0.50:
			q2 += 1
		elif j <= 0.75:
			q3 += 1
		elif j <= 1:
			q4 += 1
		else:
			continue
	return [q1, q2, q3, q4]
	
def main():
	
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data frame file containing beta values with the 1st row containing sample IDs and the 1st column containing CpG IDs.")
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
	
	printlog("Reading beta file: \"%s\"" % (options.input_file))
	data = pd.read_table(options.input_file)
	samples = data.columns[1:]

	ROUT = open(options.out_file + '.r','w')
	print ('pdf(file=\"%s\", width=10, height=10)' % (options.out_file + '.pdf'),file=ROUT)

	for s in samples:
		tmp = quarter_count(data[s])
		print ('%s <- c(%s)' % (s, ','.join([str(i) for i in tmp])), file=ROUT)
	print ("cc = rev(c('#d7191c', '#fdae61', '#a6d96a', '#1a9641'))", file=ROUT)
	print ('legend = c("beta [0.00 - 0.25]", "beta [0.25 - 0.50]", "beta [0.50 - 0.75]", "beta [0.75 - 1.00]")', file=ROUT)
	print ('nm = c(%s)' % ','.join(['"' + s + '"' for s in samples]), file=ROUT)
	print ('barplot(cbind(%s), col = cc, names.arg = nm, cex.names = 0.8, ylab = "Percentage", ylim=c(0,119), las=2, legend.text = legend)' % (','.join([s + ' * 100/sum(' + s + ')' for s in samples])), file=ROUT)
	ROUT.close()
	
	try:
		subprocess.call("Rscript " + options.out_file + '.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
		pass        	

if __name__=='__main__':
	main()	
	
