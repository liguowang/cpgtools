#!/usr/bin/env python
"""
#=========================================================================================
This program performs differential CpG analysis using beta binomial model based on 
methylation proportions (in the form of "c,n", where "c" indicates "Number of reads with
methylated C", and "n" indicates "Number of total reads". Both c and n are non-negative
integers and c <= n). Below example showing input data on 2 CpGs of 3 groups (A,B, and C)
with each group has 3 replicates:
 
cgID  A_1   A_2   A_3   B_1   B_2   B_3   C_1   C_2   C_3
CpG_1 129,170 166,178 7,9 1 6,16  10,10 10,15 11,15 16,22 20,36     
CpG_2 0,77  0,99  0,85  0,77  1,37  3,37  0,42  0,153 0,6
...

allow for covariables. 

**
you must install R package "gamlss" before running this script
https://cran.r-project.org/web/packages/gamlss/index.html
**


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
__version__="0.1.4"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

	
def main():
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="Data file containing methylation proportions (represented by \"methyl_count,total_count\", eg. \"20,30\") with the 1st row containing sample IDs (must be unique) and the 1st column containing CpG positions or probe IDs (must be unique). This file can be a regular text file or compressed file (*.gz, *.bz2) or accessible url.")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Group file defining the biological groups of each sample as well as other covariables such as gender, age.  Sample IDs should match to the \"Data file\".")
	parser.add_option("-f","--family",action="store",type="int",dest="family_func",default=1, help="A gamlss (https://cran.r-project.org/web/packages/gamlss/index.html) family object. Can be integer 1 or 2 with 1 = \"BB (beta binomial)\", 2 = \"ZIBB (zero inflated beta binomial)\" or 3 = \"ZABB (zero adjusted beta binomial)\". Default=%default.")
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
	
	FOUT = open(options.out_file + '.pval.txt','w')
	ROUT = open(options.out_file + '.r','w')
	family = {1:'BB',2:'ZIBB', 3:'ZABB'}
	if not options.family_func in family.keys():
		print ("Incorrect value of '-f'!") 
		sys.exit(106)
	
	
	print ('library("gamlss")', file=ROUT)
	
	printlog("Read group file \"%s\" ..." % (options.group_file))
	(samples,cv_names, cvs) = read_grp_file2(options.group_file)
	for cv_name in cv_names:
		print (cv_name)
		for sample in samples:
			print ('\t' + sample + '\t' + cvs[cv_name][sample])
	
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
			print ('',file=ROUT)
			print ('cgid <- \"!%s\"' % cg_id, file=ROUT)
			print ("y <- c(%s)" % (','.join([str(read) for read in methyl_reads])), file=ROUT)	#response variable
			print ("total_reads <- c(%s)" % (','.join([str(read) for read in total_reads])), file=ROUT)	#For a binomial GLM prior weights are used to give the number of trials when the response is the proportion of successes
			for cv_name in cv_names:
				print (cv_name + ' <- c(%s)' % (','.join([str(cvs[cv_name][s]) for s in  sample_IDs  ])), file = ROUT)
			
			print ('capture.output(try(fit <- gamlss(cbind(y, total_reads-y) ~ %s, family=%s)),file="NUL")' % ('+'.join(cv_names), family[options.family_func]), file = ROUT)
			print ('capture.output(s <- summary(fit, save=TRUE), file="NUL")', file=ROUT)
			print ('pval <- s$coef.table[,4]',file=ROUT)
			print ('coef <- s$coef.table[,1]',file=ROUT)
			#print ('cat(cgid, names(pval),as.vector(pval),as.vector(coef), sep="\\t")', file=ROUT)
			print ( 'write.table(file=\"%s\",x=matrix(c(cgid, names(pval),as.vector(pval),as.vector(coef)), nrow=1),quote=FALSE, row.names=FALSE, sep="\\t", col.names=FALSE, append = TRUE)' % (options.out_file + '.tmp_results.txt'),  file = ROUT) 
	ROUT.close()
	
	#sys.exit(0)
	try:
		printlog("Runing Rscript file \"%s\" ..." % (options.out_file + '.r'))
		subprocess.call("Rscript %s 2>%s" % (options.out_file + '.r', options.out_file + '.tmp_warnings.txt' ), shell=True)
	except:
		print ("Error: cannot run Rscript: \"%s\"" % (options.out_file + '.r'), file=sys.stderr)
		sys.exit(1)
	
	
	printlog("Reading file \"%s\" ..." % (options.out_file + '.tmp_results.txt'))
	bbr_results = {}
	for l in open(options.out_file + '.tmp_results.txt'):
		bbr_names = []
		bbr_pvalues = []
		bbr_coeffs = []
		l = l.strip()
		if not l.startswith('!'):continue
		l = l.replace(')','')
		l = l.replace('(','')
		f = l.split()
		cgID = f[0].replace('!','')
		tmp = f[1:]
		bbr_results[cgID] = [cv_names, ["NaN"]* len(cv_names), ["NaN"]* len(cv_names)]
		
		if len(tmp)%3 == 0:
			chunk_size = int(len(tmp)/3)
			sub_lists = [tmp[i:i+chunk_size] for i in range(0,len(tmp),chunk_size)]
			r_names = sub_lists[0]
			r_pvalues = sub_lists[1]
			r_coefs = sub_lists[2]	
			for i in range(len(r_names)):
				if r_names[i] not in cv_names:
					continue
				bbr_names.append(r_names[i])
				bbr_pvalues.append(r_pvalues[i])
				bbr_coeffs.append(r_coefs[i])
		bbr_results[cgID] = [bbr_names, bbr_pvalues, bbr_coeffs]	
	printlog("Results saved to \"%s\" ..." % (options.out_file + '.pval.txt'))
	line_num = 0
	for l in ireader.reader(options.input_file):
		line_num += 1
		f = l.split()
		if line_num == 1:
			print (l + '\t' + '\t'.join([i + '.coef' for i in bbr_names]) + '\t' + '\t'.join([i + '.pval' for i in bbr_names]), file=FOUT)
		else:
			cgID = f[0]
			print (l + '\t' + '\t'.join(bbr_results[cgID][2]) + '\t' + '\t'.join(bbr_results[cgID][1]), file=FOUT)
	
	FOUT.close()

if __name__=='__main__':
	main()
