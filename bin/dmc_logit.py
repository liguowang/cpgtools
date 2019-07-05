#!/usr/bin/env python
"""
Description
-----------
This program performs differential CpG analysis using logistic regression model based on
methylation proportions (in the form of "c,n", where "c" indicates "Number of reads with
methylated C", and "n" indicates "Number of total reads". Both c and n are  non-negative
integers and c <= n). 

Example of input data
---------------------
Below example showing input data on 2 CpGs of 3 groups (A,B, and C)
with each group has 3 replicates:
 
cgID  A_1   A_2   A_3   B_1   B_2   B_3   C_1   C_2   C_3
CpG_1 129,170 166,178 7,9 1 6,16  10,10 10,15 11,15 16,22 20,36     
CpG_2 0,77  0,99  0,85  0,77  1,37  3,37  0,42  0,153 0,6

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
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file containing methylation proportions (represented by \"methyl_count,total_count\", eg. \"20,30\") with the 1st row containing sample IDs (must be unique) and the 1st column containing CpG positions or probe IDs (must be unique). This file can be a regular text file or compressed file (*.gz, *.bz2) or accessible url..")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Group file defining the biological groups of each sample as well as other covariables such as gender, age. The first varialbe is grouping variable (must be categorical), all the other variables are considered as covariates (can be categorial or continuous). Sample IDs shoud match to the \"Data file\".")
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
	
	ROUT = open(options.out_file + '.r','w')
	family = {1:'quasibinomial', 2:'binomial',}
	if not options.family_func in family.keys():
		print ("Incorrect value of '-f'!") 
		sys.exit(106)
		
	printlog("Read group file \"%s\" ..." % (options.group_file))
	(samples,cv_names, cvs, v_types) = read_grp_file2(options.group_file)
	for cv_name in cv_names:
		print ("%s: %s" % (cv_name, v_types[cv_name]))
		for sample in samples:
			print ('\t' + sample + '\t' + cvs[cv_name][sample])
	
	print ('lrf1 <- function (cgid, m,t,%s){' % ','.join(cv_names), file=ROUT)
	print ('try(fit <- glm(cbind(m,t - m) ~ %s, family=%s))' % ('+'.join(cv_names),family[options.family_func]), file=ROUT)
	print ('pvals <- coef(summary(fit))[,4]', file=ROUT)
	print ('coefs <- coef(summary(fit))[,1]', file=ROUT)
	print ('\tif(max(pvals, na.rm=T)>1){pvals = pvals + NA}', file=ROUT)
	print ('\tif(sum(m, na.rm=T) == 0){pvals = pvals + NA}', file=ROUT)
	print ( 'write.table(file=\"%s\",x=matrix(c(cgid, as.vector(coefs), as.vector(pvals)), nrow=1),quote=FALSE, row.names=FALSE, sep="\\t", col.names=c("ID",paste(gsub("2","",names(coefs)), "coef",sep="."), paste(gsub("2","",names(pvals)), "pval",sep=".")))' % (options.out_file + '.results.txt'),  file = ROUT) 
	print ('}', file=ROUT)	
	print ('\n', file=ROUT)

	print ('lrf2 <- function (cgid, m,t,%s){' % ','.join(cv_names), file=ROUT)
	print ('try(fit <- glm(cbind(m,t - m) ~ %s, family=%s))' % ('+'.join(cv_names),family[options.family_func]), file=ROUT)
	print ('pvals <- coef(summary(fit))[,4]', file=ROUT)
	print ('coefs <- coef(summary(fit))[,1]', file=ROUT)
	print ('\tif(max(pvals, na.rm=T)>1){pvals = pvals + NA}', file=ROUT)
	print ('\tif(sum(m, na.rm=T) == 0){pvals = pvals + NA}', file=ROUT)
	print ( 'write.table(file=\"%s\",x=matrix(c(cgid, as.vector(coefs), as.vector(pvals)), nrow=1),quote=FALSE, row.names=FALSE, sep="\\t", col.names=FALSE, append = TRUE)' % (options.out_file + '.results.txt'),  file = ROUT) 
	print ('}', file=ROUT)	
	print ('\n', file=ROUT)
		
	printlog("Processing file \"%s\" ..." % (options.input_file))
	line_num = 0
	probe_list = []
	p_list = []
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
			total_reads = []	# n
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
			if line_num == 2:
				print ('lrf1(\"%s\", c(%s), c(%s), %s)' % (cg_id, ','.join([str(read) for read in methyl_reads]), ','.join([str(read) for read in total_reads]), ','.join(cv_names)), file=ROUT)
			else:
				print ('lrf2(\"%s\", c(%s), c(%s), %s)' % (cg_id, ','.join([str(read) for read in methyl_reads]), ','.join([str(read) for read in total_reads]), ','.join(cv_names)), file=ROUT)

	ROUT.close()
	
	
	try:
		printlog("Runing Rscript file \"%s\" ..." % (options.out_file + '.r'))
		subprocess.call("Rscript %s 2>%s" % (options.out_file + '.r', options.out_file + '.warnings.txt' ), shell=True)
	except:
		print ("Error: cannot run Rscript: \"%s\"" % (options.out_file + '.r'), file=sys.stderr)
		sys.exit(1)
	
if __name__=='__main__':
	main()

