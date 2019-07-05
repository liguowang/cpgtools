#!/usr/bin/env python
"""
Description
-----------
This program calculates the distribution of CpG over gene-centered genomic regions
including 'Coding exons', 'UTR exons', 'Introns', ' Upstream intergenic regions', and 
'Downsteam intergenic regions'.

Notes
-----
Please note, a particular genomic region can be assigned to different groups listed above,
because most genes have multiple transcripts, and different genes could overlap on the
genome. For example, a exon of gene A could be located in a intron of gene B. To address
this issue, we define the priority order as  below:
 0) Coding exons
 1) UTR exons
 2) Introns
 3) Upstream intergenic regions
 4) Downsteam intergenic regions  
Higher-priority group override the low-priority group. For example, if a certain part 
of a intron is overlapped with exon of other transcripts/genes, the overlapped part will 
be considered as exon (i.e. removed from intron) since "exon" has higher priority. 

#=========================================================================================
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
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="BED file specifying the C position. This BED file should have at least 3 columns (Chrom, ChromStart, ChromeEnd).  Note: the first base in a chromosome is numbered 0. This file can be a regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-r","--refgene",action="store",type="string",dest="gene_file",help="Reference gene model in standard BED-12 format (https://genome.ucsc.edu/FAQ/FAQformat.html#format1). ")
	parser.add_option("-d","--downstream",action="store",type="int",dest="downstream_size",default=2000,help="Size of down-stream intergenic region w.r.t. TES (transcription end site). default=%default (bp)")
	parser.add_option("-u","--upstream",action="store",type="int",dest="upstream_size",default=2000,help="Size of up-stream intergenic region w.r.t. TSS (transcription start site). default=%default (bp)")
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
	
	FOUT = open(options.out_file + '.tsv','w')
	ROUT = open(options.out_file + '.r','w')
	
	#step1: read CpG file
	printlog("Reading CpG file: \"%s\"" % (options.input_file))
	cpg_ranges = read_CpG_bed(options.input_file)
		
	#step2: read gene file
	printlog("Reading reference gene model: \"%s\"" % (options.gene_file))
	ref_gene = BED.ParseBED(options.gene_file)
	
	result = [("Priority_order", "Name", "Number_of_regions", "Size_of_regions(bp)", "CpG_raw_count", "CpG_count_per_KB")]
	
	#priority order: #1
	printlog("Extract Coding exons ...")
	cds_exons = ref_gene.getCDSExons(stranded=False)
	printlog("Merge Coding exons ...")
	cds_exons = BED.unionBed3(cds_exons)
	printlog("Count CpGs in Coding exons ...")
	(size,count) = count_over_range(cds_exons, cpg_ranges)
	result.append(['0','Coding exons', len(cds_exons), size, count, count*1000.0/size])	#Class, number_of_region, size_of_region, CpG_raw_count, CpG_count_perKb

	#priority order: #2
	printlog("Extract UTR exons ...")
	utr_exons = ref_gene.getUTRs(utr=35, uniquify=True, stranded = False)
	
	printlog("Merge UTR exons ...")
	utr_exons = BED.unionBed3(utr_exons)
	
	printlog("Subtract regions with higher priority from UTR exons ...")
	utr_exons = BED.subtractBed3(utr_exons, cds_exons)	#nucleotides of utr_exons that overlaps with coding exons will be removed
	
	printlog("Count CpGs in UTR exons ...")
	(size,count) = count_over_range(utr_exons, cpg_ranges)
	result.append(['1','UTR exons', len(utr_exons), size, count, count*1000.0/size])
	
	#priority order: #3
	printlog("Extract introns ...")
	introns = ref_gene.getIntrons(itype='all', uniquify=True, stranded=False)
	
	printlog("Merge introns ...")
	introns = BED.unionBed3(introns)
	
	printlog("Subtract regions with higher priority from introns ...")
	introns = BED.subtractBed3(introns, cds_exons)
	introns = BED.subtractBed3(introns, utr_exons)
	
	printlog("Count CpGs in introns ...")
	(size,count) = count_over_range(introns, cpg_ranges)
	result.append(['2','Introns', len(introns), size, count, count*1000.0/size])

	#priority order: #4
	printlog("Extract upstream intergenic regions ...")
	upstream = ref_gene.getIntergenic(direction='up', size=options.upstream_size, uniquify=True, stranded = False)
	
	printlog("Merge upstream intergenic regions ...")
	upstream = BED.unionBed3(upstream)
	
	printlog("Subtract regions with higher priority from upstream intergenic regions...")
	upstream = BED.subtractBed3(upstream, cds_exons)
	upstream = BED.subtractBed3(upstream, utr_exons)
	upstream = BED.subtractBed3(upstream, introns)
	
	printlog("Count CpGs in upstream regions ...")
	(size,count) = count_over_range(upstream, cpg_ranges)
	result.append(['3','Upstream of TSS', len(upstream), size, count, count*1000.0/size])
	
	#priority order: #5
	printlog("Extract downstream intergenic regions ...")
	downstream = ref_gene.getIntergenic(direction='down', size=options.downstream_size, uniquify=True, stranded = False)
	
	printlog("Merge downstream intergenic regions ...")
	downstream = BED.unionBed3(downstream)
	
	printlog("Subtract regions with higher priority from downstream intergenic regions...")
	downstream = BED.subtractBed3(downstream, cds_exons)
	downstream = BED.subtractBed3(downstream, utr_exons)
	downstream = BED.subtractBed3(downstream, introns)
	downstream = BED.subtractBed3(downstream, upstream)
	
	printlog("Count CpGs in downstream regions ...")
	(size,count) = count_over_range(downstream, cpg_ranges)
	result.append(['4','Downstream of TES', len(downstream), size, count, count*1000.0/size])

	print('\n')
	names=[]	#[0,1,2,3,4]
	labels = []	#[bed names]
	density=[]
	for tmp in result:
		print ('\t'.join([str(i) for i in tmp]), file=FOUT)
		names.append(tmp[0])
		labels.append(tmp[1])
		density.append(tmp[5])
	FOUT.close()
	
	print("name = c(%s)" % ','.join(['"' + i + '"' for i in names[1:]]), file=ROUT)
	print("values = c(%s)" % ','.join([str(i) for i in density[1:]]), file=ROUT)
	print ('pdf("%s", width=8, height=6)' % (options.out_file + '.pdf'), file=ROUT)
	print ('layout(matrix(c(1,1,2,1,1,2), nrow=2, byrow=TRUE))', file=ROUT)
	print ('barplot(values,names.arg=name,col=c(%s),ylab="CpG per Kb")' % ','.join(colors(5)), file=ROUT)
	print ("plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')", file=ROUT)
	for name,label in zip(names[1:], labels[1:]):
		x_pos = 0.0
		y_pos = 1-(int(name)*20.0 +5)/100 
		print ("text(x=%f, y=%f, labels=c(\"%s = %s\"),adj=c(0,0))" % (x_pos, y_pos,name,label), file=ROUT)
	print ('dev.off()', file=ROUT)
	
	ROUT.close()
	
	printlog("Running R script ...")
	try:
		subprocess.call("Rscript " + options.out_file + '.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.r', file=sys.stderr)
		pass		

if __name__=='__main__':
	main()	
	
