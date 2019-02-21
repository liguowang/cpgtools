#!/usr/bin/env python

#import built-in modules
import os,sys
import re
import string
import warnings
import string
import collections
import math
from operator import itemgetter
from itertools import groupby


#import third-party modules
from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

#from itertools import *
from cpgmodule import ireader

BED12 = '''
1. chrom - The name of the chromosome (e.g. chr3, chrY, chr2_random) or scaffold (e.g. scaffold10671).
2. chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
3. chromEnd - The ending position of the feature in the chromosome or scaffold. 
4. name - Defines the name of the BED line. 
5. score.
6. strand - Defines the strand. Either "." (=no strand) or "+" or "-".	 	 	 	 	 	 	 	 	 
7. thickStart - The starting position at which the feature is drawn thickly (for example, the start codon in gene displays). 
8. thickEnd - The ending position at which the feature is drawn thickly (for example the stop codon in gene displays).
9. itemRgb - An RGB value of the form R,G,B (e.g. 255,0,0). 
10. blockCount - The number of blocks (exons) in the BED line.
11. blockSizes - A comma-separated list of the block sizes. 
12. blockStarts - A comma-separated list of block starts.

Detailed description of BED format: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
'''


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="0.1.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"



class ParseBED:
	'''
	Manipulate BED (http://genome.ucsc.edu/FAQ/FAQformat.html) format file.
	Input BED file must be 12-column (i.e. BED-12).
	'''
	
	def __init__(self,bedFile):
		'''This is constructor of ParseBED'''
		self.f=bedFile
		self.fileName=os.path.basename(bedFile)
		self.ABS_fileName=bedFile

	def getExons(self,uniquify = True, stranded = True):
		'''
		Get all exons (including both coding exons and UTR exons) from BED-12 file.
		uniquify: if the returned blocks should be uniquify. 
		'''
		
		reblocks = []
		for l in ireader.reader(self.f):
			l = l.strip()
			if l.startswith(('#','track','browser')):continue
			f = l.split()
			if len(f) < 12:
				print ("Standard BED format has 12 columns.\n%s" % (BED), file=sys.stderr)
				sys.exit(1)
			chrom = f[0]
			chrom_start = int(f[1])
			name = f[4]
			strand = f[5]
			cdsStart = int(f[6])
			cdsEnd = int(f[7])
			blockCount = int(f[9])
			blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
			blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
			for base,offset in zip( blockStarts, blockSizes ):
				if stranded:
					reblocks.append((chrom, base, base+offset, strand))
				else:
					reblocks.append((chrom, base, base+offset))
		#self.f.seek(0)
		if uniquify:
			return list(set(reblocks))
		else:
			return reblocks
		
	def getCDSExons(self,uniquify=True, stranded = True):
		
		'''
		Get only CDS exon regions from BED-12 file. Both 5' and 3' UTR parts are removed.
		uniquify: if the returned blocks should be uniquify. 
		'''	
		reblocks = []	
		for l in ireader.reader(self.f):
			l = l.strip()
			if l.startswith(('#','track','browser')):continue
			f = l.split()
			if len(f) < 12:
				print ("\nInput error!\nStandard BED format has 12 columns.\n%s" % (BED12), file=sys.stderr)
				sys.exit(1)

			chrom = f[0]
			chrom_start = int(f[1])
			name = f[4]
			strand = f[5]
			cdsStart = int(f[6])
			cdsEnd = int(f[7])
			blockCount = int(f[9])
			blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
			blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
			cds_exons = []
			genome_seq_index = []
			for base,offset in zip( blockStarts, blockSizes ):
				if (base + offset) < cdsStart: continue
				if base > cdsEnd: continue
				exon_start = max( base, cdsStart )
				exon_end = min( base+offset, cdsEnd ) 
				if stranded:
					reblocks.append((chrom,exon_start,exon_end, strand))
				else:
					reblocks.append((chrom,exon_start,exon_end))
		#self.f.seek(0)
		if uniquify:
			return list(set(reblocks))
		else:
			return reblocks

	def getUTRs(self,utr=35, uniquify=True, stranded = True):
		'''
		Get UTR regions from BED-12 file.
		When utr=35 [default], extract both 5' and 3' UTR.
		When utr=3, only extract 3' UTR.
		When utr=5, only extract 5' UTR
		uniquify: if the returned blocks should be uniquify.
		'''
		
		reblocks = []
		for l in ireader.reader(self.f):
			l = l.strip()
			if l.startswith(('#','track','browser')):continue
			f = l.split()
			if len(f) < 12:
				print ("\nInput error!\nStandard BED format has 12 columns.\n%s" % (BED12), file=sys.stderr)
				sys.exit(1)

			chrom = f[0]
			chrom_start = int(f[1])
			name = f[4]
			strand = f[5]
			cdsStart = int(f[6])
			cdsEnd = int(f[7])
			blockCount = int(f[9])
			blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
			blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
			exon_start = []
			exon_end = []
			for base,offset in zip( blockStarts, blockSizes ):
				exon_start.append(base)
				exon_end.append(base+offset)
			if strand == '+':
				if (utr==35 or utr==5):
					for st,end in zip(exon_start,exon_end):
						if st < cdsStart:
							utr_st = st
							utr_end = min(end,cdsStart)
							if stranded:
								reblocks.append((chrom,utr_st,utr_end,strand))
							else:
								reblocks.append((chrom,utr_st,utr_end))				
				if (utr==35 or utr==3):
					for st,end in zip(exon_start,exon_end):
						if end > cdsEnd:
							utr_st = max(st, cdsEnd)
							utr_end = end
							if stranded:
								reblocks.append((chrom,utr_st,utr_end,strand))
							else:
								reblocks.append((chrom,utr_st,utr_end))				
			if strand == '-':
				if (utr==35 or utr==3):
					for st,end in zip(exon_start,exon_end):
						if st < cdsStart:
							utr_st = st
							utr_end = min(end,cdsStart)
							if stranded:
								reblocks.append((chrom,utr_st,utr_end,strand))
							else:
								reblocks.append((chrom,utr_st,utr_end))				
				if (utr==35 or utr==5):
					for st,end in zip(exon_start,exon_end):
						if end > cdsEnd:
							utr_st = max(st, cdsEnd)
							utr_end = end
							if stranded:
								reblocks.append((chrom,utr_st,utr_end,strand))
							else:
								reblocks.append((chrom,utr_st,utr_end))				
		#self.f.seek(0)
		if uniquify:
			return list(set(reblocks))
		else:
			return reblocks
				
	def getIntrons(self, itype, uniquify=True, stranded=True):
		'''
		Get Intron regions from BED-12 file. 
		separated bed file, each row represents one intron
		
		itype = :
		* 'all': all introns
		* 'first': Only return the first intron of each gene. The gene should have at least 1 intron. 
		* 'internal': return all internal introns. The gene should have at least 3 introns. 
		* 'last': Return the last intron. The gene should have at least 2 introns. 
		* 'cds': Return introns within CDS region. 
		* 'utr': Return introns within UTR regions. 
		'''

		reblocks=[]
		for l in ireader.reader(self.f):
			l = l.strip()
			if l.startswith(('#','track','browser')):continue
			f = l.split()
			chrom = f[0]
			chrom_start = int(f[1])
			name = f[4]
			strand = f[5]
			cdsStart = int(f[6])
			cdsEnd = int(f[7])
			blockCount = int(f[9])
			if blockCount == 1:continue
			blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
			blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
			exon_start = []
			exon_end = []
			for base,offset in zip( blockStarts, blockSizes ):
				exon_start.append(base)
				exon_end.append(base+offset)
	   	   
			intron_start = exon_end[:-1]
			intron_end=exon_start[1:]
			
			intron_list = list(zip(intron_start,intron_end))
			intron_number = len(intron_list)
			
			if itype == 'all':
				for (st,end) in intron_list:
					if stranded:
						reblocks.append((chrom,st,end, strand))
					else:
						reblocks.append((chrom,st,end))

			elif itype == 'first':
				if intron_number == 0:
					continue
				if strand == '-':
					if stranded:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1], strand))
					else:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1]))
				else:
					if stranded:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1], strand))
					else:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1]))
			
			elif itype == 'last':
				if intron_number < 2:
					continue
				if strand == '-':
					if stranded:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1], strand))
					else:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1]))
				else:
					if stranded:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1], strand))
					else:
						reblocks.append((chrom, intron_list[-1][0], intron_list[-1][1]))			
			elif itype == 'internal':
				if intron_number < 3:
					continue
				for (st,end) in intron_list[1:-1]:
					if stranded:
						reblocks.append((chrom,st,end, strand))
					else:
						reblocks.append((chrom,st,end))
			
			elif itype == 'cds':
				for (st,end) in  intron_list:
					if end < cdsStart: continue
					if st > cdsEnd: continue
					if stranded:
						reblocks.append((chrom,st,end, strand))
					else:
						reblocks.append((chrom,st,end))
			elif itype == 'utr':
				for (st,end) in  intron_list:
					if end < cdsStart:
						if stranded:
							reblocks.append((chrom,st,end, strand))
						else:
							reblocks.append((chrom,st,end))
					if st > cdsEnd:
						if stranded:
							reblocks.append((chrom,st,end, strand))
						else:
							reblocks.append((chrom,st,end))
									
		#self.f.seek(0)
		if uniquify:
			return list(set(reblocks))
		else:
			return reblocks


	def getIntergenic(self,direction='up', size=2000, uniquify=True, stranded = True):
		'''get intergenic regions. direction=up or down or both.'''
		
		reblocks=[]
		for l in ireader.reader(self.f):
			l = l.strip()
			if l.startswith(('#','track','browser')):continue
			f = l.split()
			chrom     = f[0]
			tx_start  = int( f[1] )
			tx_end    = int( f[2] )
			strand    = f[5]
			if(direction=="up" or direction=="both"):
				if strand=='-':
					region_st=tx_end
					region_end=tx_end +size
				else:
					region_st = max(tx_start-size,0)
					region_end=tx_start
				reblocks.append((chrom,region_st,region_end, strand))
			if (direction=="down" or direction=="both"):
				if strand == '-':
					region_st = max(0,tx_start-size)
					region_end = tx_start
				else:
					region_st = tx_end
					region_end = tx_end+size
				if stranded:
					reblocks.append((chrom,region_st,region_end, strand))
				else:
					reblocks.append((chrom,region_st,region_end))
		#self.f.seek(0)
		if uniquify:
			return list(set(reblocks))
		else:
			return reblocks

		



def unionBed3(lst):
	'''Take the union of 3 column bed files. return a new list'''
	bitsets = binned_bitsets_from_list(lst)
	ret_lst=[]
	for chrom in bitsets:
		bits = bitsets[chrom]
		end = 0
		while 1:
			start = bits.next_set( end )
			if start == bits.size: break
			end = bits.next_clear( start )
			ret_lst.append([chrom, start, end])
	bitsets=dict()
	return ret_lst

def intersectBed3(lst1,lst2):
	'''Take the intersection of two bed files (3 column bed files)'''
	bits1 = binned_bitsets_from_list(lst1)
	bits2 = binned_bitsets_from_list(lst2)

	bitsets = dict()
	ret_lst = []
	for key in bits1:
		if key in bits2:
			bits1[key].iand( bits2[key] )
			bitsets[key] = bits1[key]

	for chrom in bitsets:
		bits = bitsets[chrom]
		end = 0
		while 1:
			start = bits.next_set( end )
			if start == bits.size: break
			end = bits.next_clear( start )
			ret_lst.append([chrom, start, end])
	bits1.clear()
	bits2.clear()
	bitsets.clear()
	return ret_lst

def subtractBed3(lst1,lst2):
	'''subtrack lst2 from lst1'''
	bitsets1 = binned_bitsets_from_list(lst1)
	bitsets2 = binned_bitsets_from_list(lst2)
	
	ret_lst=[]
	for chrom in bitsets1:  
		if chrom not in bitsets1:
			continue
		bits1 = bitsets1[chrom]
		if chrom in bitsets2:
			bits2 = bitsets2[chrom]
			bits2.invert()
			bits1.iand( bits2 )
		end=0
		while 1:
			start = bits1.next_set( end )
			if start == bits1.size: break
			end = bits1.next_clear( start )
			ret_lst.append([chrom,start,end])
	bitsets1 = dict()
	bitsets2 = dict()
	return ret_lst

def tillingBed(chrName,chrSize,stepSize=10000):
	'''tilling whome genome into small sizes'''
	#tilling genome
	for start in xrange(0,chrSize,stepSize):
		end = start + stepSize
		if end < chrSize:
			yield (chrName,start,end)
		else:
			yield (chrName,start,chrSize)
		
