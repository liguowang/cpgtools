#!/usr/bin/env python
"""
Description
-----------
This program uses Bayesian Gaussian Mixture model (BGMM) to trichotmize beta values into 
three status: 
 * Un-methylated (labeled as "0" in result file)
 * Semi-methylated (labeled as "1" in result file)
 * Full-methylated (labeled as "2" in result file)
 * unassigned (labeled as "-1" in result file)
"""


import sys,os
import collections
import numpy as np
from optparse import OptionParser
from sklearn import mixture
from time import strftime
from cpgmodule import ireader
from cpgmodule.utils import *
import pandas as pd

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def load_data(infile):
	"""
	Input file is tab or space separated plain text file.
	*The first row contains sample IDs (must be unique)
	*The first column contains probe IDs (must be unique)
	*Each cell (except for the 1st row and 1st column) contains Beta-value
	
	Example:
	
	Probe			sample_1	sample_2	sample_3 ...
	cg09835024		0.0547		0.1187		0.0625	...
	cg25813447		0.428		0.3746		0.0666	...
	cg07779434		0.3713		0.4194		0.0493	...
	"""

	printlog("Reading input file: \"%s\"" % infile)
	df1 = pd.read_table(infile, index_col = 0)

	#remove any rows with NAs
	df2 = df1.dropna(axis=0, how='any')
	printlog("%d rows with missing values were removed." % (len(df1) - len(df2)))
        
	print ("\tTotal samples: %d" % (len(df2.columns)), file=sys.stderr)
	print ("\tTotal probes: %d" % len(df2), file=sys.stderr)
	return df2

def build_GMM(d,rnd):
	"""
	Return means of components of Gaussian Mixture Model.
	d is data frame returned by "load_data" function.
	rnd is a random number. You get exactly the same results when running multiple times using the same random number. Must be integer. 
	"""
	
	bgmm_models = collections.defaultdict(list)
	for s_id in sorted(d.columns):
		printlog ("Building Bayesian Gaussian Mixture model for subject: %s ...\r" % s_id)
		bgmm = mixture.BayesianGaussianMixture(n_components=3, covariance_type='full',max_iter=50000,tol=0.001,random_state=rnd)
		bgmm_models[s_id] = bgmm.fit(d[s_id].values.reshape(-1,1))
	#print (bgmm_models)
	return bgmm_models
	

def summary_GMM(m):
	"""
	Summarize BGMM models returned by "build_GMM"
	"""
	printlog ("Summerzie GMM models ...")
	
	FOUT = open('summary_report.txt','w')
	
	print ("\n\n#means of components", file=FOUT)
	print ("Subject_ID\tUnmethyl\tSemiMethyl\tMethyl",file=FOUT)
	for k,v in m.items():
		print (k + '\t' + '\t'.join([str(i) for i in sorted(v.means_[:,0])]),file=FOUT)		
	
	
	print ("\n\n#Weights of components", file=FOUT)
	print ("Subject_ID\tUnmethyl\tSemiMethyl\tMethyl", file=FOUT)
	for k,v in m.items():
		print (k + '\t' + '\t'.join([str(i) for i in sorted(v.weights_)]), file=FOUT)		
	
	
	print ("\n\n#Converge status and n_iter", file=FOUT)
	print ("Subject_ID\tConverged\tn_iter", file=FOUT)
	for k,v in m.items():
		print (k + '\t' + '\t'.join([str(i) for i in (v.converged_, v.n_iter_)]), file=FOUT)		
	FOUT.close()
	
	printlog ("Reports were saved into \"summary_report.txt\".")
	
def trichotmize(d,m, prob_cutoff):
	"""
	trichotmize beta-value into one of ('0','0.5','1')
	0 : Un-methylation
	0.5: Semi-methylation
	1: Methylation
	
	d is beta value object returned by "load_data" function
	m is BGMM models returned by 'build_GMM' function
	
	"""
	probe_IDs = list(d.index)
	
	for s_id in sorted(m.keys()):
		printlog ("Writing to \"%s\" ..." % (s_id + ".results.txt"))
		FOUT = open(s_id + ".results.txt",'w')
		methyl_lables = {}	#key is index (index can be 0,1 or 2 corresponding to 3 components), value is 0, 1 or 2 corresponding to Un-, Semi- and full- methylation
		component_means = m[s_id].means_[:,0]	# list of component means
		betas = d[s_id]
		for idx,val in enumerate(component_means):
			if val == max(component_means):
				methyl_lables[idx] = '2'	# full methyl
			elif val == min(component_means):
				methyl_lables[idx] = '0'	# un-methyl
			else:
				methyl_lables[idx] = '1'	# semi-methyl
				
		probs = m[s_id].predict_proba(d[s_id].values.reshape(-1,1))	# list of probabilities of components: [[  4.33638063e-035   9.54842259e-001   4.51577411e-002],...]
		
		print ("#Prob_of_0: Probability of CpG belonging to un-methylation group", file=FOUT)
		print ("#Prob_of_1: Probability of CpG belonging to semi-methylation group", file=FOUT)
		print ("#Prob_of_2: Probability of CpG belonging to full-methylation group", file=FOUT)
		print ("#Assigned_lable: -1 = 'unassigned', 0 = 'un-methylation', 1 = 'semi-methylation', 2 = 'full-methylation'", file=FOUT)
		print ("Probe_ID" + '\tBeta_value\t' + '\t'.join(['Prob_of_' + methyl_lables[0], 'Prob_of_' + methyl_lables[1], 'Prob_of_' + methyl_lables[2]]) + '\t' + 'Assigned_lable', file=FOUT)
		for probe_ID, beta, p in zip(probe_IDs, betas, probs):
			p_list = list(p)
			#print probe_ID
			#print p_list
			if methyl_lables[p_list.index(max(p_list))] == '1':
				if max(p_list) >= prob_cutoff:
					print (probe_ID + '\t' + str(beta) + '\t' + '\t'.join([str(i) for i in p_list]) + '\t' + methyl_lables[p_list.index(max(p_list))], file=FOUT)
				else:
					print (probe_ID + '\t' + str(beta) + '\t' + '\t'.join([str(i) for i in p_list]) + '\t' + '-1', file=FOUT)
			else:
				print (probe_ID + '\t' + str(beta) + '\t' + '\t'.join([str(i) for i in p_list]) + '\t' + methyl_lables[p_list.index(max(p_list))], file=FOUT)
		FOUT.close()

def main():
	print (__doc__)
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Input plain text file containing beta values with the 1st row containing sample IDs (must be unique) and the 1st column containing probe IDs (must be unique).")
	parser.add_option("-c","--prob-cut",action="store",type="float",dest="prob_cutoff",default=0.99,help="Probability cutoff to assign a probe into \"semi-methylated\" class. default=%default")
	parser.add_option("-r","--report",action="store_true",dest="report_summary",help="Presense of this flag renders program to generate \"summary_report.txt\" file.")
	parser.add_option("-s","--seed",action="store",type='int', dest="random_state",default=99, help="The seed used by the random number generator. default=%default")
	(options,args)=parser.parse_args()
	
	print ()

	if not (options.input_file):
		parser.print_help()
		sys.exit(0)
	if not os.path.exists(options.input_file):
		print ('\n\n' + options.input_file + " does NOT exists" + '\n',file=sys.stderr)
		sys.exit(0)
	
	infile = options.input_file
	
	#step1: read beta value file
	dat = load_data(infile)	
	
	#step2: build BGMM models
	GMMs = build_GMM(dat, rnd = options.random_state)
	
	#step3: Summerize BGMM models
	if options.report_summary:
		summary_GMM(GMMs)
	
	#step4: Classification
	trichotmize(dat, GMMs, options.prob_cutoff)
	

if __name__=='__main__':
	main()	
	
