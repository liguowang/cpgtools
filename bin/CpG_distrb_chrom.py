#!/usr/bin/env python
"""
Description
-----------
This program calculates the distribution of CpG over chromosomes.
"""


import sys,os
import collections
import subprocess
import numpy as np
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *

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
	parser.add_option("-i","--input-files",action="store",type="string",dest="input_files",help="Input CpG file(s) in BED3+ format. Multiple BED files should be separated by \",\" (eg: \"-i file_1.bed,file_2.bed,file_3.bed\"). BED file can be a regular text file or compressed file (.gz, .bz2). The barplot figures will NOT be generated if you provide more than 12 samples (bed files). [required]")
	parser.add_option("-n","--names",action="store",type="string",dest="file_names",help="Shorter and meaningful names to label samples. Should be separated by \",\" and match CpG BED files in number. If not provided, basenames of CpG BED files will be used to label samples. [optional]")
	parser.add_option("-s","--chrom-size",action="store",type="string",dest="chrom_size",help="Chromosome size file. Tab or space separated text file with 2 columns: the first column is chromosome name/ID, the second column is chromosome size. This file will determine: (1) which chromosomes are included in the final barplots, so do NOT include 'unplaced', 'alternative' contigs in this file. (2) The order of chromosomes in the final barplots.  [required]")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file. [required]")
	(options,args)=parser.parse_args()
	
	print ()

	if not (options.input_files):
		print (__doc__)
		#print ('You must specify input file(s)',file=sys.stderr)
		parser.print_help()
		sys.exit(101)
	if not (options.chrom_size):
		print (__doc__)
		#print ('You must specify the chrom size file',file=sys.stderr)
		parser.print_help()
		sys.exit(102)
	if not (options.out_file):
		print (__doc__)
		#print ('You must specify the output file',file=sys.stderr)
		parser.print_help()
		sys.exit(103)	
		
	input_files = 	options.input_files.split(',')
	for i in input_files + [options.chrom_size]:
		if not os.path.exists(i):
			print ('\n' + i + " does NOT exists" + '\n',file=sys.stderr)
			sys.exit(104)
	
	input_names = []
	if options.file_names:
		input_names = options.file_names.split(',')
	else:
		for f in input_files:
			input_names.append(os.path.basename(f))
	if len(input_files) != len(input_names):
		print ('-i and -n don\'t match in number',file=sys.stderr)
		sys.exit(105)

		
	#step1: read chrom sizes files
	printlog("Reading chromosome size file: \"%s\"" % (options.chrom_size))
	cnames,csizes = read_chromSize(options.chrom_size)
	for cname,csize in zip(cnames,csizes):
		print("  " + cname + '\t' + str(csize))
	
	#step2: read CpG files
	dat = collections.defaultdict(dict)	#samleName:chromName:CpGount
	for f,n in zip(input_files, input_names):
		printlog("Reading CpG BED file \"%s\" named \"%s\"" % (f,n))
		dat[n] = chrom_count(f)
	
	
	#step3: write matrix to file
	printlog("Save CpG count to \"%s\"" % (options.out_file + '.txt'))
	FOUT = open(options.out_file + '.txt','w')
	print ("chromID\tchromSize\t" + '\t'.join([ n + '.CpG_count' for n in input_names]), file=FOUT)
	
	for cname,csize in zip(cnames,csizes):
		row = []
		row.append(cname.replace('chr',''))
		row.append(csize)
		for n in input_names:
			try:
				row.append(dat[n][cname])
			except:
				row.append(0)
		print ('\t'.join([str(i) for i in row]), file=FOUT)
		
	FOUT.close()
		
	
	#step 4: print R script
	if len(input_names) <= 12:
		printlog("Generate R script, save to \"%s\"" % (options.out_file + '.r'))
		ROUT = open(options.out_file + '.r','w')
		print ("chromNames = c(%s)" % (','.join(['"' + i.replace('chr','') + '"' for i in cnames])),file=ROUT)
		print ("chromSizes = c(%s)" % (','.join([str(i) for i in csizes])),file=ROUT)
	
		input_names2 = ['X_' + i for i in input_names]
		for n1,n2 in zip(input_names, input_names2):
			tmp = []
			for cname in cnames:
				try:
					tmp.append(dat[n1][cname])
				except:
					tmp.append(0)
			print ("%s = c(%s)" % (n2, ','.join([str(i) for i in tmp])), file=ROUT)
	
		my_col = colors(len(input_names))
		print ('cols = c(%s)' % ','.join(my_col),file=ROUT)
	
	
		print ('pdf(\"%s\", width=12, height=6)' %  (options.out_file + '.CpG_total.pdf'), file=ROUT)
		print ('barplot(rbind(%s),col=cols,beside=T,names.arg=%s, xlab="Chromosome", ylab="CpG count", legend.text=c(%s), cex.names=0.5, cex.axis=0.6)' % (','.join(input_names2), 'chromNames', ','.join(['"' + i + '"' for i in input_names])), file=ROUT)	
		print ('dev.off()', file=ROUT)
	
		print ('pdf(\"%s\", width=12, height=6)' %  (options.out_file + '.CpG_percent.pdf'), file=ROUT)
		print ('barplot(rbind(%s),col=cols,beside=T,names.arg=%s, xlab="Chromosome", ylab="CpG percent", legend.text=c(%s), cex.names=0.5, cex.axis=0.6)' % (','.join([ i + '*100.0/sum(' + i + ')' for i in input_names2]), 'chromNames', ','.join(['"' + i + '"' for i in input_names])), file=ROUT)	
		print ('dev.off()', file=ROUT)
	
		print ('pdf(\"%s\", width=12, height=6)' %  (options.out_file + '.CpG_perMb.pdf'), file=ROUT)
		print ('barplot(rbind(%s),col=cols,beside=T,names.arg=%s, xlab="Chromosome", ylab="CpG per Mb", legend.text=c(%s), cex.names=0.5, cex.axis=0.6)' % (','.join([ i + '*1000000.0/chromSizes' for i in input_names2]), 'chromNames', ','.join(['"' + i + '"' for i in input_names])), file=ROUT)	
		print ('dev.off()', file=ROUT)
	
		ROUT.close()

		#step 5: Run R script
		printlog("Running R script ...")
		try:
			subprocess.call("Rscript " + options.out_file + '.r', shell=True)
		except:
			print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
			pass
	else:
		print ("Cannot generate R script file and pdf files.", file=sys.stderr)
    	
    	            


if __name__=='__main__':
	main()	
	
