#!/usr/bin/env python
"""
#=========================================================================================
This program performs differential CpG analysis using logistic regression model based on
methylation proportions (in the form of "c,n", where "c" indicates "Number of reads with
methylated C", and "n" indicates "Number of total reads". Both c and n are  non-negative
integers and c <= n). Below example showing input data on 2 CpGs of 3 groups (A,B, and C)
with each group has 3 replicates:
 
cgID  A_1   A_2   A_3   B_1   B_2   B_3   C_1   C_2   C_3
CpG_1 129,170 166,178 7,9 1 6,16  10,10 10,15 11,15 16,22 20,36     
CpG_2 0,77  0,99  0,85  0,77  1,37  3,37  0,42  0,153 0,6

allow for covariables. 
...

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
__version__="0.1.6"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

	
def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file containing methylation proportions (represented by \"methyl_count,total_count\", eg. \"20,30\") with the 1st row containing sample IDs (must be unique) and the 1st column containing CpG positions or probe IDs (must be unique). This file can be a regular text file or compressed file (*.gz, *.bz2) or accessible url..")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Group file defining the biological groups of each sample as well as other covariables such as gender, age. The first varialbe is usually categorical and used to make the contrast (calculate pvalues), all the other variables are considered as covariates.   Sample IDs shoud match to the \"Data file\".")
	parser.add_option("-f","--family",action="store",type="int",dest="family_func",default=1, help="Error distribution and link function to be used in the GLM model. Can be integer 1 or 2 with 1 = \"quasibinomial\" and 2 = \"binomial\". Default=%default.")
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
	if not os.path.isfile(options.input_file):
		print ("Input data file \"%s\" does not exist\n" % options.input_file) 
		sys.exit(104)
	if not os.path.isfile(options.group_file):
		print ("Input group file \"%s\" does not exist\n" % options.input_file) 
		sys.exit(105)
	if os.path.exists(options.out_file + '.results.txt'):
		os.remove(options.out_file + '.results.txt')
		
	ROUT = open(options.out_file + '.r','w')
	family = {1:'quasibinomial',2:'binomial'}
	if not options.family_func in family.keys():
		print ("Incorrect value of '-f'!") 
		sys.exit(106)
		
	printlog("Read group file \"%s\" ..." % (options.group_file))
	(samples,cv_names, cvs, v_types) = read_grp_file2(options.group_file)
	for cv_name in cv_names:
		print ("%s: %s" % (cv_name, v_types[cv_name]))
		for sample in samples:
			print ('\t' + sample + '\t' + cvs[cv_name][sample])
	
	print ('lrf <- function (cgid, m,t,%s){' % ','.join(cv_names), file=ROUT)
	print ('\ttry(fit1 <- glm(cbind(m,t - m) ~ %s, family=%s))' % ('+'.join(cv_names),family[options.family_func]), file=ROUT)
	if len(cv_names) == 1:
		print ('\ttry(fit0 <- glm(cbind(m,t - m) ~ 1, family=%s))' % (family[options.family_func]), file=ROUT)
	elif len(cv_names) >1:
		print ('\ttry(fit0 <- glm(cbind(m,t - m) ~ %s, family=%s))' % ('+'.join(cv_names[1:]),family[options.family_func]), file=ROUT)
	print ('\ttest <- anova(fit1, fit0,test="Chisq")', file=ROUT)
	print ('\tpval <- test$"Pr(>Chi)"[2]', file=ROUT)
	print ('\tresults <- list(cgID = cgid, pvalue = pval)', file=ROUT)
	print ('\twrite.table(file=\"%s\",x=results, quote=FALSE, row.names=FALSE, sep="\\t",append = TRUE, col.names=FALSE)' % (options.out_file + '.results.txt'),  file = ROUT)
	print ('}', file=ROUT)	
	print ('\n', file=ROUT)

		
	printlog("Processing file \"%s\" ..." % (options.input_file))
	line_num = 0
	for l in ireader.reader(options.input_file):
		line_num += 1
		f = l.split()
		if line_num == 1:
			sample_IDs = f[1:]
			# check if sample ID matches
			for s in samples:
				if s not in sample_IDs:
					printlog("Cannot find sample ID \"%s\" from file \"%s\"" % (s, options.input_file))
					sys.exit(3)
			#for cv_name in cv_names:
			#	print (cv_name + ' <- c(%s)' % (','.join([str(cvs[cv_name][s]) for s in  sample_IDs  ])), file = ROUT)
			for cv_name in cv_names:
				if v_types[cv_name] == 'continuous':
					print (cv_name + ' <- c(%s)' % (','.join([str(cvs[cv_name][s]) for s in  sample_IDs  ])), file = ROUT)
				elif  v_types[cv_name] == 'categorical':
					print (cv_name + ' <- as.factor(c(%s))' % (','.join([str(cvs[cv_name][s]) for s in  sample_IDs  ])), file = ROUT)
				else:
					printlog("unknown vaiable type!")
					sys.exit(1)

			print ('\n', file=ROUT)
			continue
		else:
			methyl_reads = []			# c
			total_reads = []			# n
			cg_id = f[0]
			for i in f[1:]:
				#try:
				m = re.match(r'(\d+)\s*\,\s*(\d+)', i)
				if m is None:
					methyl_reads.append("NaN")
					total_reads.append("NaN")
					continue
				else:
					c = int(m.group(1))
					n = int(m.group(2))
					if n >= c and n > 0:
						methyl_reads.append(c)
						total_reads.append(n)
					else:
						printlog("Incorrect data format!")
						print (f)
						sys.exit(1)		
			print ('lrf(\"%s\", c(%s), c(%s), %s)' % (cg_id, ','.join([str(read) for read in methyl_reads]), ','.join([str(read) for read in total_reads]), ','.join(cv_names)), file=ROUT)

	ROUT.close()
	
	
	try:
		printlog("Runing Rscript file \"%s\" ..." % (options.out_file + '.r'))
		subprocess.call("Rscript %s 2>%s" % (options.out_file + '.r', options.out_file + '.r.warnings.txt' ), shell=True)
	except:
		print ("Error: cannot run Rscript: \"%s\"" % (options.out_file + '.r'), file=sys.stderr)
		sys.exit(1)
	
	
	printlog("Perfrom Benjamini-Hochberg (aka FDR) correction ...")
	probe_list0 = []	#probes without valid pvalue
	probe_list1 = []
	p_list1 = []	
	if os.path.exists(options.out_file + '.results.txt') and os.path.getsize(options.out_file + '.results.txt') > 0:
		for l in ireader.reader(options.out_file + '.results.txt'):
			f = l.split()
			id = f[0]
			try:
				pv = float(f[1])
				probe_list1.append(id)
				p_list1.append(pv)
			except:
				probe_list0.append(id)
				continue
	q_list1 =  padjust.multiple_testing_correction(p_list1)
	
	OUT = open(options.out_file + '.results.txt','w')
	print ("probe\tP-value\tadj.Pvalue", file = OUT)
	
	#probes with valid p and q
	for id,p,q in zip(probe_list1, p_list1, q_list1):
		print (id + '\t' + str(p) + '\t' + str(q), file=OUT)
	
	#probes without valid p and q
	if len(probe_list0) > 0:
		for id in probe_list0:
			print (id + '\tNA\tNA', file=OUT)
	
	
	OUT.close()
	
if __name__=='__main__':
	main()

