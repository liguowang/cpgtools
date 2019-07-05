#!/usr/bin/env python
"""
Description
-----------
This program calculates the methylation profile (i.e. average beta value) for genomic regions 
around genes. These genomic regions include: "5'UTR exon", "CDS exon", "3'UTR exon",
"first intron", "internal intron", "last intron", "up-stream intergenic", and
"down-stream intergenic".

Example of input BED6+ file
---------------------------
chr22   44021512        44021513        cg24055475      0.9231  -
chr13   111568382       111568383       cg06540715      0.1071  +
chr20   44033594        44033595        cg21482942      0.6122  -
"""


import sys,os
import collections
import subprocess
import numpy as np
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED

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
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="BED6+ file specifying the C position. This BED file should have at least 6 columns (Chrom, ChromStart, ChromeEnd, Name, Beta_value, Strand). BED6+ file can be a regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-r","--refgene",action="store",type="string",dest="gene_file",help="Reference gene model in standard BED12 format (https://genome.ucsc.edu/FAQ/FAQformat.html#format1). \"Strand\" column must exist in order to decide 5' and 3' UTRs, up- and down-stream intergenic regions.")
	parser.add_option("-d","--downstream",action="store",type="int",dest="downstream_size",default=2000,help="Size of down-stream genomic region added to gene. default=%default (bp)")
	parser.add_option("-u","--upstream",action="store",type="int",dest="upstream_size",default=2000,help="Size of up-stream genomic region added to gene. default=%default (bp)")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file.")
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
	
	FOUT = open(options.out_file + '.txt','w')
	ROUT = open(options.out_file + '.r','w')
	print ("\t".join(["Group","Relative_position(5'->3')", "Average_beta"]), file=FOUT)
	
	#step1: read CpG file
	printlog("Reading CpG file: \"%s\"" % (options.input_file))
	cpg_ranges = read_CpG_bed(options.input_file)
		
	#step2: read gene file
	printlog("Reading reference gene model: \"%s\"" % (options.gene_file))
	ref_gene = BED.ParseBED(options.gene_file)
	
	group_sizes = []	#number of datapoints in each group
	printlog("Process upstream regions ...")
	up_2k = ref_gene.getIntergenic(direction = 'up',  size=options.upstream_size)
	s = coverage_over_range(up_2k,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['Upstream_intergenic',str(i), str(s[i])]), file=FOUT)	
	print ('Upstream_intergenic_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)
	

	printlog("Process 5' UTR exons ...")
	utr5_exons = ref_gene.getUTRs(utr=5)
	s = coverage_over_range(utr5_exons,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['Five_prime_UTR',str(i), str(s[i])]), file=FOUT)	
	print ('Five_prime_UTR_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)

	
	printlog("Process Coding exons ...")
	cds_exons = ref_gene.getCDSExons()
	s = coverage_over_range(cds_exons,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['Coding_exon',str(i), str(s[i])]), file=FOUT)
	print ('Coding_exon_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)

	printlog("Process first introns ...")
	introns = ref_gene.getIntrons(itype='first')
	s = coverage_over_range(introns,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['First_intron',str(i), str(s[i])]), file=FOUT)
	print ('First_intron_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)
			
	printlog("Process internal introns ...")
	introns = ref_gene.getIntrons(itype='internal')
	s = coverage_over_range(introns,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['Internal_intron',str(i), str(s[i])]), file=FOUT)
	print ('Internal_intron_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)


	printlog("Process last introns ...")
	introns = ref_gene.getIntrons(itype='last')
	s = coverage_over_range(introns,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['Last_intron',str(i), str(s[i])]), file=FOUT)
	print ('Last_intron_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)
	
					
	printlog("Process 3' UTR exons ...")
	utr3_exons = ref_gene.getUTRs(utr=3)
	s = coverage_over_range(utr3_exons,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['Three_prime_UTR',str(i), str(s[i])]), file=FOUT)
	print ('Three_prime_UTR_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)

	printlog("Process downstream regions ...")
	down_2k = ref_gene.getIntergenic(direction = 'down', size=options.downstream_size)
	s = coverage_over_range(down_2k,cpg_ranges)
	group_sizes.append(len(s))
	for i in sorted(s):
		print ('\t'.join(['Downstream_intergenic',str(i), str(s[i])]), file=FOUT)
	print ('Downstream_intergenic_y=c(%s)' % ','.join([str(s[i]) for i in sorted(s)]), file=ROUT)
	
	print('\n')
	print ('pdf(file=\"%s\", width=10, height=5)' % (options.out_file + '.pdf'),file=ROUT)
	print ('plot(1:%d, c(Upstream_intergenic_y, Five_prime_UTR_y, Coding_exon_y, First_intron_y, Internal_intron_y, Last_intron_y, Three_prime_UTR_y, Downstream_intergenic_y),ylim=c(0,1), xaxt="n",xlab="", ylab="Average methylation", type="l", col="red")' % sum(group_sizes), file=ROUT)
	print ('abline(v = c(100,201,302,403,504,605,706),col="blue", lty="dashed")', file=ROUT)
	print ('abline(v = c(%d,%d,%d,%d,%d,%d,%d),col="blue", lty="dashed")' % (sum(group_sizes[0:1]),sum(group_sizes[0:2]),sum(group_sizes[0:3]),sum(group_sizes[0:4]),sum(group_sizes[0:5]),sum(group_sizes[0:6]),sum(group_sizes[0:7])), file=ROUT)
	print ('abline(h = 0.5,col="grey", lty="dashed")', file=ROUT)
	print ('text(x=c(%d,%d,%d,%d,%d,%d,%d, %d)+50, y=0.9, cex=0.7, labels=c("Upstream\\n(5\'->3\')", "5\'UTR exon\\n(5\'->3\')","Coding exon\\n(5\'->3\')","First intron\\n(5\'->3\')","Internal intron\\n(5\'->3\')","Last intron\\n(5\'->3\')", "3\'UTR exon\\n(5\'->3\')","Downstream\n(5\'->3\')"))' % (0, sum(group_sizes[0:1]),sum(group_sizes[0:2]),sum(group_sizes[0:3]),sum(group_sizes[0:4]),sum(group_sizes[0:5]),sum(group_sizes[0:6]),sum(group_sizes[0:7])), file=ROUT)
	print ('dev.off()',file=ROUT)
		
	FOUT.close()
	ROUT.close()
	try:
		subprocess.call("Rscript " + options.out_file + '.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
		pass        


if __name__=='__main__':
	main()	
	
