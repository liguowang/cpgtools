#!/usr/bin/env python
"""
Description
-----------
This program performs differential CpG analysis using Fisher's exact test. It only applies
to two sample comparison with no replicates. if replicates are provided, *methyl reads*
and *total reads* of all replicates will be sumed

Example of input data file
--------------------------
cgID        sample_1    sample_2
CpG_1       129,170     166,178
CpG_2       24,77       67,99 

number before "," indicates *number of methyl reads*
number after "," indicates *number of total reads*

Output
-------
Three columns ("Odds ratio", "pvalue" and "adjusted pvalue") will append to input data table. 
#=========================================================================================
"""


import sys,os
import collections
import subprocess
import numpy as np
import re
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

	
def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file containing methylation proportions (represented by \"methyl_count,total_count\", eg. \"20,30\") with the 1st row containing sample IDs (must be unique) and the 1st column containing CpG positions or probe IDs (must be unique). This file can be a regular text file or compressed file (*.gz, *.bz2) or accessible url.")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Group file defining the biological groups of each sample. It is a comma-separated 2 columns file with the 1st column containing sample IDs, and the 2nd column containing group IDs.  It must have a header row. Sample IDs should match to the \"Data file\".")
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
	(s,g) = read_grp_file1(options.group_file)
	s2g = dict(zip(s,g))
	g2s = collections.defaultdict(list)
	
	for k,v in s2g.items():
		g2s[v].append(k)
	
	group_IDs = sorted(g2s.keys())
	for g in group_IDs:
		print ("\tGroup %s has %d samples:" % (g, len(g2s[g])))
		print ('\t\t' + ','.join(g2s[g]))
	
	if len(group_IDs) != 2:
		printlog("You must have two groups!", file=sys.stderr)
		sys.exit(1)
	
	line_num = 1
	probe_list = []
	p_list = []
	or_list = []
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
			cg_id = f[0]
			probe_list.append(cg_id)
			proportions = f[1:]
			methyl_reads = 0
			unmethyl_reads = 0
			g2values = collections.defaultdict(dict)
			for g in group_IDs:
				g2values[g]['methyl'] = 0
				g2values[g]['unmethyl'] = 0
			for s,p in zip(sample_IDs, proportions):
				gid = s2g[s]
				m = re.match(r'(\d+)\s*\,\s*(\d+)', p)
				if m is None:
					continue
				else:
					c = int(m.group(1))
					n = int(m.group(2))
					if n >= c and n > 0:
						g2values[gid]['methyl'] += c
						g2values[gid]['unmethyl']  += (n-c)
					else:
						printlog("Incorrect data format!")
						print (f)
						sys.exit(1)		
			(odds, pval) = stats.fisher_exact([ [g2values[group_IDs[0]]['methyl'], g2values[group_IDs[0]]['unmethyl']],[g2values[group_IDs[1]]['methyl'], g2values[group_IDs[1]]['unmethyl']] ])
			#print (g2values[group_IDs[0]]['methyl'], g2values[group_IDs[0]]['unmethyl'],g2values[group_IDs[1]]['methyl'], g2values[group_IDs[1]]['unmethyl'])
			p_list.append(pval)
			or_list.append(odds)				
		line_num += 1
	
	printlog("Perfrom Benjamini-Hochberg (aka FDR) correction ...")
	adjusted_p = {}
	q_list =  padjust.multiple_testing_correction(p_list)
	for id,o,p,q in zip(probe_list, or_list, p_list, q_list):
		adjusted_p[id] = '\t'.join([str(i) for i in (o,p,q)])
	
	printlog("Writing to %s" % (options.out_file + '.pval.txt'))
	line_num = 1
	for l in ireader.reader(options.input_file):
		if line_num == 1:
			print (l + '\tOddsRatio\tpval\tadj.pval', file=FOUT)
		else:
			f = l.split()
			probe_ID = f[0]
			print (l + '\t' + adjusted_p[probe_ID], file=FOUT)
		line_num += 1
	FOUT.close()
		
	
	

if __name__=='__main__':
	main()
