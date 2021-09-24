import sys
from bx.intervals import *
import numpy as np
from cpgmodule import ireader

def getBasalDomains(bedfile, printit = False):
	'''
	Define gene's basal regulatory domain. 
	bedfile: one gene one TSS (could use the canonical (longest) isoform, or merge all isoforms into a super transcript.
	'''
	basal_ranges = {}
	
	for l in ireader.reader(bedfile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		try:
			chrom = f[0]
			start = int(f[1])
			end = int(f[2])
			symbol = f[3]
			gene_strand = f[5]
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue
			if gene_strand not in ['+','-']:
				print ("Invalid strand. Skip: " + l, file=sys.stderr)
				continue					
		except:
			print ("BED has at lesat 6 columns. Skip: " + l, file=sys.stderr)
			continue

		if chrom not in basal_ranges:
			basal_ranges[chrom] = IntervalTree()
		
		basal_ranges[chrom].insert_interval( Interval(start, end, strand = gene_strand, value = symbol))
		
		if printit:
			print('\t'.join([str(i) for i in (chrom, basal_st, basal_end, symbol, '0', strand)]), file = sys.stdout)
	return basal_ranges

def geteExtendedDomains(basal_ranges, bedfile, up_ext=2000, down_ext=2000, min_gene = 200, printit = False):
	'''
	Define gene's extended regulatory domain. 
	bedfile:one gene one TSS (could use the canonical (longest) isoform, or merge all
			isoforms into a super transcript.
	up_ext:
		Size of extension to upstream. Should be multiples of 100
	down_ext:
		Size of extension to downstream. Should be multiples of 100
	min_gene:
		minimum gene size (from TSS to TES). Should be multiples of 100
	
	'''	
	return_ranges = []
	
	for l in ireader.reader(bedfile):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		try:
			chrom = f[0]
			start = int(f[1])
			end = int(f[2])
			symbol = f[3]
			strand = f[5]
			
			if start < 0:continue
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue	
			if (end - start ) < min_gene:
				continue		
			if strand not in ['+', '-']:
				print ("Unknown strand. Skip: " + l, file=sys.stderr)
				continue
		except:
			print ("BED has at lesat 6 columns. Skip: " + l, file=sys.stderr)

		if strand == '+':
			extension_st = start  - up_ext
			extension_end = end + down_ext
		elif strand == '-':
			extension_st = start  - down_ext
			extension_end = end + up_ext
		if extension_st < 0: 
			extension_st = 0
		
		#try to update extension_st
		overlaps = basal_ranges[chrom].find(extension_st, start)
		if len(overlaps) > 0:
			for o in overlaps:
				if o.end > extension_st:
					extension_st = o.end
			if extension_st > start:
				extension_st = start
		
		if (start - extension_st) < min_gene:
			continue
		
		#try to update extension_end
		overlaps = basal_ranges[chrom].find(end, extension_end)
		if len(overlaps) > 0:
			for o in overlaps:
				if o.start < extension_end:
					extension_end = o.start
			if extension_end < end:
				extension_end = end
		
		if (extension_end - end) < min_gene:
			continue
		
		return_ranges.append(([chrom, extension_st, start,symbol], [chrom, start, end,symbol], [chrom, end, extension_end,symbol], strand))
		#return_ranges.append(([chrom, extension_st, start, symbol], [chrom, start, end, symbol], [chrom, end, extension_end,symbol], strand))
		#return_ranges.append(([chrom, extension_st, start, strand], [chrom, start, end, strand], [chrom, end, extension_end, strand]))

		if printit:
			print('\t'.join([str(i) for i in (chrom, extension_st, extension_end, symbol, '0', strand,  start, end, '255,0,0', 1, extension_end - extension_st, 0)]), file = sys.stdout)
		
	return return_ranges

if __name__=='__main__': 		
	tmp = getBasalDomains(sys.argv[1], printit = False)	
	b = geteExtendedDomains(basal_ranges = tmp, bedfile = sys.argv[1], printit=False)	
	
	
	for a1,a2,a3,a4 in b:
		if a4 == '+':
			print ('\t'.join([str(i) for i in a1]) + '(UIR)\t' +str(int(a1[2]) - int(a1[1])) +  '\t+')
			print ('\t'.join([str(i) for i in a2]) + '(Body)\t' +str(int(a2[2]) - int(a2[1])) +  '\t+')
			print ('\t'.join([str(i) for i in a3]) + '(DIR)\t' +str(int(a3[2]) - int(a3[1])) +  '\t+')
		if a4 == '-':
			print ('\t'.join([str(i) for i in a1]) + '(DIR)\t' +str(int(a1[2]) - int(a1[1])) +  '\t-')
			print ('\t'.join([str(i) for i in a2]) + '(Body)\t' +str(int(a2[2]) - int(a2[1])) +  '\t-')
			print ('\t'.join([str(i) for i in a3]) + '(UIR)\t' +str(int(a3[2]) - int(a3[1])) +  '\t-')
		
		
		
		
