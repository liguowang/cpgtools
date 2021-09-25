#!/usr/bin/env python3

"""
#=========================================================================================
Select the K best features according to the K highest scores. Scores can be measured by:

* ANOVA F-value between label/feature for classification tasks.
* Mutual information for a discrete target.	
* Chi-squared stats of non-negative features for classification tasks.

Example of input data file
---------------------------
CpG_ID	Sample_01	Sample_02	Sample_03	Sample_04
cg_001	0.831035	0.878022	0.794427	0.880911
cg_002	0.249544	0.209949	0.234294	0.236680
cg_003	0.845065	0.843957	0.840184	0.824286
"""
import sys
import numpy as np
from optparse import OptionParser
from cpgmodule.utils import *
import pandas as pd

from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import chi2,f_classif,mutual_info_classif

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.10.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

	
def main():
	
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input_file",action="store",type="string",dest="input_file",help="Tab-separated data frame file containing beta values with the 1st row containing sample IDs and the 1st column containing CpG IDs.")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Comma-separated group file defining the biological groups of each sample.")
	parser.add_option("-k","--topK",action="store",type='int', dest="cpg_count", default=100, help="Number of top features to select. default=%default" )
	parser.add_option("-s","--score-function",action="store",type='string', dest="score_function", default='chi2', help="Scoring function used to measure the dependency between features scores and labels. Must be \"chisq\" (chi-squared statistic), \"anova\" (ANOVA F-value), or \"mutual_info\" (mutual information). default=%default" )
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="The prefix of the output file.")
	(options,args)=parser.parse_args()
	
	print ()
	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)
	if not (options.group_file):
		print (__doc__)
		parser.print_help()
		sys.exit(102)	
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)	
	
	printlog("Reading input file: \"%s\"" % (options.input_file))
	df1 = pd.read_csv(options.input_file, index_col = 0, sep="\t")
	#print (df1)
	
	#remove any rows with NAs
	df2 = df1.dropna(axis=0, how='any')
	printlog("%d rows with missing values were removed." % (len(df1) - len(df2)))
	#print (df2)
	
	printlog("Transposing data matrix ... ")
	df2 = df2.T
	total_feature = len(df2.columns)
	printlog("Total number of features: %d " % (total_feature))
	#print (df2)
	
	
	printlog("Reading group file: \"%s\"" % (options.group_file))
	group = pd.read_csv(options.group_file, index_col=0, header=0,names=['Sample_ID', 'Group_ID'])
	a = list(group['Group_ID'])	#a is *string labels* for groups: ['Normal', 'Normal', 'Normal', 'Normal', 'Normal', 'Cancer', 'Cancer', 'Cancer', 'Cancer']
	#print (a)
	y, tmp = pd.factorize(a)	#y is *numeric labels* for groups: [0 0 0 0 0 1 1 1 1]
	#print (np.array(y))
	
	if options.cpg_count < total_feature:
		
		if options.score_function == 'anova':
			printlog ("Using ANOVA F value to select features ...")
			selector = SelectKBest(f_classif, k = options.cpg_count)
		elif options.score_function == 'mutual_info':
			printlog ("Using Mutual Information to select features ...")
			selector = SelectKBest(mutual_info_classif, k = options.cpg_count)
		else:
			printlog ("Using Chi Square statistic to select features ...")
			selector = SelectKBest(chi2, k = options.cpg_count)
	else:
		printlog("Doing nothing! '-k' >= the total number of features in \"%s\"" % (options.input_file))
		sys.exit(0)
	
	
	selector.fit_transform(df2, np.array(y))
	cols = selector.get_support(indices=False)
	selected_data = df2.loc[:,cols]
	selected_featureNum = len(selected_data.columns)
	printlog("Total number of selected features : %d " % (selected_featureNum))
	#print (selected_data)
	
	printlog("Writing to file: \"%s\"" % (options.out_file + '.selectedFeatures.tsv'))
	pd.DataFrame.to_csv(selected_data.T, options.out_file + '.selectedFeatures.tsv', sep="\t", index_label="sample")
	

if __name__=='__main__':
	main()	