# CpGtools
# Copyright (c) 2024-2026 Liguo Wang
#
# Author: Liguo Wang
# Email: "deep.omics.lab@gmail.com"
# Project: https://github.com/liguowang/cpgtools
#
# This file is part of CpGtools and is distributed under the MIT License.
# See the LICENSE.txt file in the project root for the full license text.

import os,sys

from bx.bitset import *
from bx.bitset_builders import *
from bx.intervals import *

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
			#name = f[4]
			strand = f[5]
			cdsStart = int(f[6])
			cdsEnd = int(f[7])
			#blockCount = int(f[9])
			blockSizes = [ int(i) for i in f[10].strip(',').split(',') ]
			blockStarts = [ chrom_start + int(i) for i in f[11].strip(',').split(',') ]
			#cds_exons = []
			#genome_seq_index = []
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
			#name = f[4]
			strand = f[5]
			cdsStart = int(f[6])
			cdsEnd = int(f[7])
			#blockCount = int(f[9])
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
	    """
	    Extract intron regions from a BED12 transcript annotation.

	    Parameters
	    ----------
	    itype : {"all", "first", "internal", "last", "cds", "utr"}
	        Type of introns to return.

	        * ``all``:
	          Return all introns.
	        * ``first``:
	          Return the first intron in transcript orientation.
	          Transcripts must contain at least one intron.
	        * ``internal``:
	          Return introns excluding the first and last introns.
	          Transcripts must contain at least three introns.
	        * ``last``:
	          Return the last intron in transcript orientation.
	          Transcripts must contain at least two introns.
	        * ``cds``:
	          Return introns overlapping the coding region.
	        * ``utr``:
	          Return introns located entirely outside the coding region.

	    uniquify : bool, optional
	        Remove duplicate intervals before returning them.

	    stranded : bool, optional
	        Include strand as the fourth field of each returned interval.

	    Returns
	    -------
	    list
	        Intron intervals as ``(chrom, start, end)`` or
	        ``(chrom, start, end, strand)`` tuples.

	    Notes
	    -----
	    BED12 exon blocks are stored in ascending genomic-coordinate order.
	    For minus-strand transcripts, the intron list is reversed before selecting
	    first, internal, and last introns so that these categories are defined in
	    transcript 5'-to-3' orientation.
	    """
	    valid_types = {"all", "first", "internal", "last", "cds", "utr"}
	    if itype not in valid_types:
	        raise ValueError(
	            f"Unsupported intron type: {itype!r}. "
	            f"Choose from: {', '.join(sorted(valid_types))}."
	        )

	    reblocks = []

	    def append_interval(chrom, start, end, strand):
	        if start >= end:
	            return

	        if stranded:
	            reblocks.append((chrom, start, end, strand))
	        else:
	            reblocks.append((chrom, start, end))

	    for line in ireader.reader(self.f):
	        line = line.strip()

	        if not line or line.startswith(("#", "track", "browser")):
	            continue

	        fields = line.split()
	        if len(fields) < 12:
	            raise ValueError(
	                "Input gene annotation must be in BED12 format. "
	                f"Found {len(fields)} columns in line: {line}"
	            )

	        chrom = fields[0]
	        chrom_start = int(fields[1])
	        strand = fields[5]
	        cds_start = int(fields[6])
	        cds_end = int(fields[7])
	        block_count = int(fields[9])

	        if strand not in {"+", "-"}:
	            raise ValueError(
	                f"Invalid transcript strand {strand!r}; expected '+' or '-'."
	            )

	        if block_count < 2:
	            continue

	        block_sizes = [
	            int(value)
	            for value in fields[10].rstrip(",").split(",")
	            if value
	        ]
	        relative_block_starts = [
	            int(value)
	            for value in fields[11].rstrip(",").split(",")
	            if value
	        ]

	        if (
	            len(block_sizes) != block_count
	            or len(relative_block_starts) != block_count
	        ):
	            raise ValueError(
	                "BED12 blockCount does not match blockSizes/blockStarts "
	                f"for transcript {fields[3]!r}."
	            )

	        exon_starts = [
	            chrom_start + relative_start
	            for relative_start in relative_block_starts
	        ]
	        exon_ends = [
	            start + size
	            for start, size in zip(exon_starts, block_sizes)
	        ]

	        genomic_introns = [
	            (left_exon_end, right_exon_start)
	            for left_exon_end, right_exon_start
	            in zip(exon_ends[:-1], exon_starts[1:])
	            if left_exon_end < right_exon_start
	        ]

	        if not genomic_introns:
	            continue

	        # Convert genomic order to transcript 5'-to-3' order.
	        transcript_introns = (
	            genomic_introns
	            if strand == "+"
	            else list(reversed(genomic_introns))
	        )

	        if itype == "all":
	            selected_introns = genomic_introns

	        elif itype == "first":
	            selected_introns = [transcript_introns[0]]

	        elif itype == "last":
	            if len(transcript_introns) < 2:
	                continue
	            selected_introns = [transcript_introns[-1]]

	        elif itype == "internal":
	            if len(transcript_introns) < 3:
	                continue
	            selected_introns = transcript_introns[1:-1]

	        elif itype == "cds":
	            selected_introns = [
	                (start, end)
	                for start, end in genomic_introns
	                if start < cds_end and end > cds_start
	            ]

	        else:  # itype == "utr"
	            selected_introns = [
	                (start, end)
	                for start, end in genomic_introns
	                if end <= cds_start or start >= cds_end
	            ]

	        for start, end in selected_introns:
	            append_interval(chrom, start, end, strand)

	    if uniquify:
	        return list(dict.fromkeys(reblocks))

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
		
