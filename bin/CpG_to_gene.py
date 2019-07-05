#!/usr/bin/env python
"""
Description
-----------
This program annotates CpGs by assigning them to their putative target genes. Follows the 
"Basel plus extension" rules used by GREAT(http://great.stanford.edu/public/html/index.php)
 
 * Basal regulatory domain:
   is a user-defined genomic region around the TSS (transcription start site). By default,
   from TSS upstream 5kb to TSS downstream 1Kb is considered as the gene's *basal regulatory
   domain*. When defining a gene's "basal regulatory domain", the other nearby genes will be
   ignored (which means different genes' basal regulatory domains can be overlapped.)

 * Extended regulatory domain:
   The gene regulatory domain is extended in both directions to the nearest gene's "basal 
   regulatory domain" but no more than the maximum extension (default = 1000 kb) in one
   direction.

Notes
-----
 1. Which genes are assigned to a particular CpG largely depends on gene annotation. A
    "conservative" gene model (such as Refseq curated protein coding genes) is recommended. 
 2. In the gene model, multiple isoforms should be merged into a single gene.
#=========================================================================================
"""


import sys,os
import collections
import subprocess
import numpy as np
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule.region2gene import *

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
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="BED3+ file specifying the C position. BED3+ file could be a regular text file or compressed file (.gz, .bz2). [required]")
	parser.add_option("-r","--refgene",action="store",type="string",dest="gene_file",help="Reference gene model in BED12 format (https://genome.ucsc.edu/FAQ/FAQformat.html#format1). \"One gene one transcript\" is recommended. Since most genes have multiple transcripts; one can collapse multiple transcripts of the same gene into a single super transcript or select the canonical transcript.")
	parser.add_option("-u","--basal-up",action="store",type="int",dest="basal_up_size",default=5000,help="Size of extension to upstream of TSS (used to define gene's \"basal regulatory domain\"). default=%default (bp)")
	parser.add_option("-d","--basal-down",action="store",type="int",dest="basal_down_size",default=1000,help="Size of extension to downstream of TSS (used to define gene's basal regulatory domain). default=%default (bp)")
	parser.add_option("-e","--extension",action="store",type="int",dest="extension_size",default=1000000,help="Size of extension to both up- and down-stream of TSS (used to define gene's \"extended regulatory domain\"). default=%default (bp)")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of the output file. Two additional columns will be appended to the original BED file with the last column indicating \"genes whose extended regulatory domain are overlapped with the CpG\", the 2nd last column indicating \"genes whose basal regulatory domain are overlapped with the CpG\". [required]")
	(options,args)=parser.parse_args()
	
	print ()

	if not (options.input_file):
		#print ('You must specify input file(s)',file=sys.stderr)
		print (__doc__)
		parser.print_help()
		sys.exit(101)
	if not (options.gene_file):
		#print ('You must specify the chrom size file',file=sys.stderr)
		print (__doc__)
		parser.print_help()
		sys.exit(102)
	if not (options.out_file):
		#print ('You must specify the output file',file=sys.stderr)
		print (__doc__)
		parser.print_help()
		sys.exit(103)	
	
	FOUT = open(options.out_file + '.associated_genes.txt','w')
	print ("#The last column contains genes whose extended regulatory domain are overlapped with the CpG", file=FOUT)
	print ("#The 2nd last column contains genes whose basal regulatory domain are overlapped with the CpG", file=FOUT)
	print ("#\"//\" indicates no genes are found", file=FOUT)
	
	printlog("Calculate basal regulatory domain from: \"%s\" ..." % (options.gene_file))
	basal_domains = getBasalDomains(bedfile = options.gene_file, up = options.basal_up_size, down = options.basal_down_size, printit = False)
	
	printlog("Calculate extended regulatory domain from: \"%s\" ..." % (options.gene_file))
	extended_domains = geteExtendedDomains(basal_ranges = basal_domains, bedfile = options.gene_file, up = options.basal_up_size, down = options.basal_down_size, ext=options.extension_size, printit = False)
	
	#overlap = extended_domains['chr1'].find(2161048,2161049)
	
	printlog("Assigning CpG to gene ...")
	for l in ireader.reader(options.input_file):
		if l.startswith('#'):
			print (l, file=FOUT)
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		try:
			f = l.split()
			chrom = f[0]	
			start = int(f[1])
			end = int(f[2])
		except:
			print ("Invalid BED line: %s" % l, file=sys.stderr)
			continue
						
		
		basal_genes = set()	#genes whose basal domain is overlapped with CpG
		if chrom not in basal_domains:
			basal_genes.add('//')
		else:
			overlaps = basal_domains[chrom].find(start,end)
			if len(overlaps) == 0:
				basal_genes.add('//')
			else:
				for o in overlaps:
					basal_genes.add(o.value)
		
		extend_genes = set()	#genes whose extended domain is overlapped with CpG
		if chrom not in extended_domains:
			extend_genes.add('//')
		else:
			overlaps = extended_domains[chrom].find(start,end)
			if len(overlaps) == 0:
				extend_genes.add('//')
			else:
				for o in overlaps:
					extend_genes.add(o.value)
		
		
		extend_genes = extend_genes - basal_genes
		if len(extend_genes) == 0:
			extend_genes.add('//')
		print (l + '\t' + ';'.join(basal_genes) + '\t' + ';'.join(extend_genes), file=FOUT)
	FOUT.close()

if __name__=='__main__': 
	main()	
		
	
