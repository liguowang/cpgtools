#!/usr/bin/env python3

"""
Description
-----------
This program calculates the CpG density (count) profile over gene body as well as its up-
down-stream regions. It is useful to visualize how CpGs are distributed around genes.

Specifically, the up-stream region, gene region (from TSS to TES) and down-stream region
wil be equally divided into 100 bins, then CpG count was aggregated over a total of 300 bins
from 5' to 3' (upstream bins, gene bins, downstrem bins).
#==========================================================================================
"""

import sys
import subprocess
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED
from cpgmodule import extend_bed

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
	parser.add_option("-i","--input_file",action="store",type="string",dest="input_file",help="BED file specifying the C position. This BED file should have at least three columns (Chrom, ChromStart, ChromeEnd).  Note: the first base in a chromosome is numbered 0. This file can be a regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-r","--refgene",action="store",type="string",dest="gene_file",help="Reference gene model in standard BED6+ format.")
	parser.add_option("-d","--downstream",action="store",type="int",dest="downstream_size",default=2000,help="Maximum extension size from TES (transcription end site) to down-stream to define the \"downstream intergenic region (DIR)\". Note: (1) The actual used DIR size can be smaller because the extending process could stop earlier if it reaches the boundary of another nearby gene. (2) If the actual used DIR size is smaller than cutoff defined by \"-c/--SizeCut\", the gene will be skipped.  default=%default (bp)")
	parser.add_option("-u","--upstream",action="store",type="int",dest="upstream_size",default=2000,help="Maximum extension size from TSS (transcription start site) to up-stream to define the \"upstream intergenic region (UIR)\". Note: (1) The actual used UIR size can be smaller because the extending process could stop earlier if it reaches the boundary of another nearby gene. (2) If the actual used UIR size is smaller than cutoff defined by \"-c/--SizeCut\", the gene will be skipped. default=%default (bp)")
	parser.add_option("-c","--SizeCut",action="store",type="int",dest="minimum_size",default=200,help="The minimum gene size. Gene size is defined as the genomic size between TSS and TES, including both exons and introns. default=%default (bp)")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="The prefix of the output file.")
	(options,args)=parser.parse_args()

	print ()

	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)

	if not (options.gene_file):
		print (__doc__)
		parser.print_help()
		sys.exit(102)

	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)

	FOUT = open(options.out_file + '.tsv','w')
	ROUT = open(options.out_file + '.r','w')

	#step1: read CpG file
	printlog("Reading CpG file: \"%s\"" % (options.input_file))
	cpg_ranges = read_CpG_bed(options.input_file)

	#step2: read gene file
	printlog("Reading reference gene model: \"%s\"" % (options.gene_file))
	tmp1 = extend_bed.getBasalDomains(options.gene_file)
	tmp2 = extend_bed.geteExtendedDomains(basal_ranges = tmp1, bedfile = options.gene_file, up_ext=options.upstream_size, down_ext=options.downstream_size, min_gene = options.minimum_size, printit = False)

	printlog("Calculating CpG density ...")
	#CpG density
	(up_density, gene_density, down_density) = density_over_range(tmp2, cpg_ranges)

	printlog("Wrting data to : \"%s\"" % (options.out_file + '.tsv'))
	print ("Group\tPosition\tCpG_count", file=FOUT)
	for ind in (sorted(up_density)):
		print ("Upstream\t" + str(ind) + '\t' + str(up_density[ind]), file = FOUT)

	for ind in (sorted(gene_density)):
		print ("GeneBody\t" + str(ind) + '\t' + str(gene_density[ind]), file = FOUT)

	for ind in (sorted(down_density)):
		print ("Downstream\t" + str(ind) + '\t' + str(down_density[ind]), file = FOUT)
	FOUT.close()

	print ('pdf(file=\"%s\", width=10, height=5)' % (options.out_file + '.pdf'),file=ROUT)
	print ("d <- read.table(file = '%s', header = T, sep='\\t')" % (options.out_file + '.tsv'), file = ROUT)
	print ("x = 1:length(d$CpG_count)", file=ROUT)
	print ("plot(x,d$CpG_count,type='l',col='red',lwd=1,xaxt='n',ylab='CpG count',xlab='')", file=ROUT)
	print ("abline(v = c(102,203),col='blue', lty='dashed', lwd=0.5)", file=ROUT)
	print ("text(x=c(0,102,203)+50, y=0.1, labels=c('Upstream', 'geneBody','Downstream'))", file=ROUT)
	print ('dev.off()',file=ROUT)
	ROUT.close()

	printlog("Running R script to: '%s'" % (options.out_file + '.r'))
	try:
		subprocess.call("Rscript " + options.out_file + '.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
		pass

if __name__=='__main__':
	main()

