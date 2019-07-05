#!/usr/bin/env python
"""
Description
-----------
This program performs differential CpG analysis based on beta values. It uses Student's
t-test for two group comparison, and ANOVA for multiple groups comparison.
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

def standard_ttest(a, b, equalVar=True, nanPolicy='omit'):
	'''
	Calculate the T-test for the means of two independent samples of scores.
	'''
	p = np.nan
	t = np.nan
	try:
		tmp = stats.ttest_ind(a, b, equal_var = equalVar, nan_policy = nanPolicy)
		p = tmp.pvalue
		t = tmp.statistic
	except:
		pass
	return (p,t)

def paired_ttest(a, b, nanPolicy='omit'):
	'''
	Calculate the T-test on TWO RELATED samples of scores, a and b.
	'''
	p = np.nan
	t = np.nan
	try:
		tmp = stats.ttest_rel(a,b, nan_policy = nanPolicy)
		p = tmp.pvalue
		t = tmp.statistic
	except:
		pass
			
	return (p,t)

def anova(*args):
	'''
	The one-way ANOVA tests the null hypothesis that three or more groups have the same population mean
	'''
	p = np.nan
	t = np.nan
	try:
		tmp = stats.f_oneway(*args)
		p = tmp.pvalue
		t = tmp.statistic
	except:
		pass
	return (p,t)
	
def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file containing beta values with the 1st row containing sample IDs (must be unique) and the 1st column containing CpG positions or probe IDs (must be unique). Except for the 1st row and 1st column, any non-numerical values will be considered as \"missing values\" and ignored. This file can be a regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Group file defining the biological group of each sample. It is a comma-separated 2 columns file with the 1st column containing sample IDs, and the 2nd column containing group IDs.  It must have a header row. Sample IDs should match to the \"Data file\". Note: automatically switch to use ANOVA if more than 2 groups were defined in this file.")
	parser.add_option("-p","--paired",action="store_true",default=False,dest="paired",help="If '-p/--paired' flag was specified, use paired t-test which requires the equal number of samples in both groups. Paired sampels are matched by the order. This option will be ignored for multiple group analysis.")
	parser.add_option("-w","--welch",action="store_true",default=False,dest="welch_ttest",help="If '-w/--welch' flag was specified, using Welch's t-test which does not assume the two samples have equal variance.  If omitted, use standard two-sample t-test (i.e. assuming the two samples have equal variance). This option will be ignored for paired t-test and multiple group analysis.")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
	(options,args)=parser.parse_args()
	
	print ()
	#print (options.paired)
	#print (options.welch_ttest)
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
	#ROUT = open(options.out_file + '.r','w')
	
	printlog("Read group file \"%s\" ..." % (options.group_file))
	(ss,gs) = read_grp_file1(options.group_file)
	
	s2g = {}
	for s,g in zip(ss,gs):
		s2g[s] = g	
	
	g2s = collections.defaultdict(list)
	for s,g in zip(ss, gs):
		g2s[g].append(s)
	
	group_IDs = sorted(g2s.keys())
	for g in group_IDs:
		print ("\tGroup %s has %d samples:" % (g, len(g2s[g])))
		print ('\t\t' + ','.join(g2s[g]))
	
	if len(group_IDs) < 2:
		printlog("You must have at least two groups!", file=sys.stderr)
		sys.exit(1)
	elif (len(group_IDs) == 2) and (options.paired is True):
		printlog("Perfrom paired t-test of two related samples ...")
		if len(g2s[group_IDs[0]]) != len(g2s[group_IDs[1]]):
			printlog("Unequal sample size. Cannot perform paired t-test.")
			sys.exit(2)
	elif (len(group_IDs) == 2) and (options.paired is False):
		printlog("Perfrom standard t-test of two independent samples ...")
	elif len(group_IDs) >= 3:
		printlog("Perfrom ANOVA ...")
	
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
				if options.paired:
					(pval,tscore) = paired_ttest(a,b)
				else:
					(pval,tscore) = standard_ttest(a,b, equalVar = options.welch_ttest)				
			elif len(g2values) >= 3:
				tmp = []
				for g in group_IDs:
					tmp.append(np.array(g2values[g]))
				(pval,tscore) = anova(*tmp)
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
