import sys
from bx.intervals import *
import numpy as np
from cpgmodule import ireader

def getBasalDomains(bedfile, up = 5000, down = 1000, printit = False):
	'''
	Define gene's basal regulatory domain. 
	bedfile: one gene one TSS (could use the canonical (longest) isoform, or merge all isoforms into a super transcript.
	up: size of extension to upstream of TSS
	down: size of extension to downstream of TSS
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
			strand = f[5]
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue		
		except:
			print ("BED has at lesat 6 columns. Skip: " + l, file=sys.stderr)
			continue

		if chrom not in basal_ranges:
			basal_ranges[chrom] = IntervalTree()
		
		if strand == '+':
			basal_st = (start + 1) - up
			basal_end = (start + 1) + down
			basal_st = max(0, basal_st)
			basal_ranges[chrom].insert_interval( Interval(basal_st, basal_end, strand = '+', value = symbol))

		elif strand == '-':
			basal_st = end - down
			basal_end = end + up
			basal_st = max(0, basal_st)
			basal_ranges[chrom].insert_interval( Interval(basal_st, basal_end, strand = '-', value = symbol))
		if printit:
			print('\t'.join([str(i) for i in (chrom, basal_st, basal_end, symbol, '0', strand)]), file = sys.stdout)
	return basal_ranges

def geteExtendedDomains(basal_ranges, bedfile, up = 5000, down = 1000, ext=1000000, printit = False):
	'''
	Define gene's extended regulatory domain. 
	bedfile:one gene one TSS (could use the canonical (longest) isoform, or merge all
			isoforms into a super transcript.
	ext: 	maximum size of extension (default 1000Kb)
	
	Two step process:
	1) Each gene is assigned a basal regulatory domain of a minimum distance upstream and
	   downstream of the TSS (regardless of other nearby genes). 
	2) The gene regulatory domain is extended in both directions to the nearest gene's
	    basal domain but no more than the maximum extension in one direction.
	'''	
	domain_ranges = {}	#gene's regulatory domain range
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
			if start > end:
				print ("'Start' cannot be larger than 'End'. Skip: " + l, file=sys.stderr)
				continue		
		except:
			print ("BED has at lesat 6 columns. Skip: " + l, file=sys.stderr)

		if strand == '+':
			tss = start + 1
			basal_st = tss - up
			basal_end = tss + down
			extension_st = tss - ext
			extension_end = tss + ext
		elif strand == '-':
			tss = end
			basal_st = tss - down
			basal_end = tss + up
			extension_st = tss - ext
			extension_end = tss + ext
		if extension_st < 0: 
			extension_st = 0
		#try to update extension_st
		overlaps = basal_ranges[chrom].find(extension_st, basal_st)
		if len(overlaps) > 0:
			for o in overlaps:
				if o.end > extension_st:
					extension_st = o.end
			if extension_st > basal_st:
				extension_st = basal_st

		#try to update extension_end
		overlaps = basal_ranges[chrom].find(basal_end, extension_end)
		if len(overlaps) > 0:
			for o in overlaps:
				if o.start < extension_end:
					extension_end = o.start
			if extension_end < basal_end:
				extension_end = basal_end
								
		if chrom not in domain_ranges:
			domain_ranges[chrom] = IntervalTree()
		else:
			domain_ranges[chrom].insert_interval(Interval(extension_st, extension_end, strand = strand, value=symbol))


		if printit:
			print('\t'.join([str(i) for i in (chrom, extension_st, extension_end, symbol, '0', strand,  basal_st, basal_end, '255,0,0', 1, extension_end - extension_st, 0)]), file = sys.stdout)
		
	return domain_ranges
	"""
		if len(overlaps) == 1:
			domain_ranges[chrom].insert_interval(Interval(extension_st, extension_end, strand = strand, value=symbol))
			if printit:
				print('\t'.join([str(i) for i in (chrom, extension_st, extension_end, symbol, '0', strand)]), file = sys.stdout)

		else:
			o_basal_starts = []	#starts of basal_region overlapped with extension_region
			o_basal_ends = []
			for o in overlaps:
				if o.vallue = symbol:
					continue
				o_basal_starts.append(o.start)
				o_basal_ends.append(o.end)
			
			tmp1 = [i for i in o_basal_ends if i > extension_st and i < tss]
			tmp2 = [i for i in o_basal_starts if i < extension_end and i > tss]
			if len(tmp1) == 0:
				truncaed_ext_st = extension_st
			else:
				truncaed_ext_st  = max(tmp1)
			if len(tmp2) == 0:
				truncaed_ext_end = extension_end
			else:
				truncaed_ext_end = min(tmp2)
			
			truncaed_ext_st = max(0,truncaed_ext_st)
			domain_ranges[chrom].insert_interval(Interval(truncaed_ext_st, truncaed_ext_end, strand = strand, value=symbol))
			
			if printit:
				print('\t'.join([str(i) for i in (chrom, truncaed_ext_st, truncaed_ext_end, symbol + '_extended', '0', strand)]), file = sys.stdout)
	"""
								
if __name__=='__main__': 		
	tmp = basal_domain(sys.argv[1], printit = False)	
	extended_domain(basal_ranges = tmp, bedfile = sys.argv[1], printit=True)	
		
		
		
		
		
