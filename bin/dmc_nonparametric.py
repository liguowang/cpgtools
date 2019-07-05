#!/usr/bin/env python
"""
Description
-----------
This program performs differential CpG analysis uisng the  Mann-Whitney U test
for two group comparison, and the Kruskal-Wallis H-test for multiple groups
comparison. 
"""


import sys,os
import collections
import subprocess
import numpy as np
from scipy import stats
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED
from cpgmodule import padjust

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def mwu_test(a, b):
	'''
	mann-whitney U test of two samples. 
	'''
	p = np.nan
	t = np.nan
	try:
		tmp = stats.mannwhitneyu(a, b,  alternative='two-sided')
		p = tmp.pvalue
		t = tmp.statistic
	except:
		pass
	return (p,t)

def kruskal_test(*args):
	'''
	Compute the Kruskal-Wallis H-test for independent samples
	'''
	p = np.nan
	t = np.nan
	try:
		tmp = stats.kruskal(*args, nan_policy='omit')
		p = tmp.pvalue
		t = tmp.statistic
	except:
		pass
	return (p,t)
	
def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file containing beta values with the 1st row containing sample IDs (must be unique) and the 1st column containing CpG positions or probe IDs (must be unique). Except for the 1st row and 1st column, any non-numerical values will be considered as \"missing values\" and ignored. This file can be a regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Group file defining the biological group of each sample. It is a comma-separated 2 columns file with the 1st column containing sample IDs, and the 2nd column containing group IDs. It must have a header row. Sample IDs should match to the \"Data file\". Note: automatically switch to use  Kruskal-Wallis H-test if more than 2 groups were defined in this file.")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
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
	
	FOUT = open(options.out_file + '.pval.txt','w')
	
	printlog("Read group file \"%s\" ..." % (options.group_file))
	(s,g) = read_grp_file1(options.group_file)
	s2g = dict(zip(s,g))
	g2s = collections.defaultdict(list)
	
	for k,v in s2g.items():
		g2s[v].append(k)
	
	group_IDs = sorted(g2s.keys())
	for g in group_IDs:
		print ("\tGroup %s has %d samples:" % (g, len(g2s[g])))
		print ('\t\t' + ','.join(g2s[g]))
	
	if len(group_IDs) < 2:
		printlog("You must have at least two groups!", file=sys.stderr)
		sys.exit(1)
	elif len(group_IDs) == 2:
		printlog("Perfrom Mann-Whitney rank test of two samples ...")
	elif len(group_IDs) >= 3:
		printlog("Perfrom Kruskal-Wallis H-test ...")
	
	line_num = 1
	probe_list = []
	p_list = []
	for l in ireader.reader(options.input_file):
		f = l.split()
		if line_num == 1:
			sample_IDs = f[1:]

			# check if sample ID matches
			for s in s2g:
				if s not in sample_IDs:
					printlog("Cannot find sample ID \"%s\" from file \"%s\"" % (s, options.input_file))
					sys.exit(3)
		else:
			g2values = collections.defaultdict(list)
			probe_ID = f[0]
			beta_values = f[1:]
			for s,b in zip(sample_IDs, beta_values):
			
				#deal with non-numerical values
				try:
					b = float(b)
				except:
					b = np.nan
				
				#skip if s not in group file
				if s not in s2g:
					continue
				
				gid = s2g[s]
				g2values[gid].append(b)			
			
			if len(g2values) == 2:
				a = np.array(g2values[group_IDs[0]])
				b = np.array(g2values[group_IDs[1]])
				(pval,tscore) = mwu_test(a,b)
			elif len(g2values) >= 3:
				tmp = []
				for g in group_IDs:
					tmp.append(np.array(g2values[g]))
				(pval,tscore) = kruskal_test(*tmp)
			probe_list.append(probe_ID)
			p_list.append(pval)
		line_num += 1
	
	printlog("Perfrom Benjamini-Hochberg (aka FDR) correction ...")
	adjusted_p = {}
	q_list =  padjust.multiple_testing_correction(p_list)
	for id,p,q in zip(probe_list, p_list, q_list):
		adjusted_p[id] = '\t'.join([str(i) for i in (p,q)])
	
	printlog("Writing to %s" % (options.out_file + '.pval.txt'))
	line_num = 1
	for l in ireader.reader(options.input_file):
		if line_num == 1:
			print (l + '\tpval\tadj.pval', file=FOUT)
		else:
			f = l.split()
			probe_ID = f[0]
			print (l + '\t' + adjusted_p[probe_ID], file=FOUT)
		line_num += 1
	FOUT.close()
		
	
	

if __name__=='__main__':
	main()
