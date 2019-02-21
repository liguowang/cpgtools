#!/usr/bin/env python
'''DNA/protein motif visualization and scan'''

#import built-in modules
import sys,os
from collections import defaultdict
from scipy import stats 
import itertools

#import third-party modules
import numpy as np
#changes to the paths

#changing history to this module


__author__ = "Liguo Wang"
__copyright__ = ""
__credits__ = []
__license__ = "GPLv2"
__version__ = "1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "Wang.Liguo@mayo.edu"
__status__ = "Development" #Prototype or Production

class PSSM (object):
	'''
	Description: provides functions to manipulate Position-Specific Scoring Matrix (PSSM) 
	such as PFM, PPM, PWM matrix.
	'''
	
	def __init__(self, sites, dna = True, name = None, rv = False):
		'''
		Initialize object.
		Must be DNA or protein (dna = False) sequences.
		dna = True: DNA sequence
		dna = False: protein sequence
		rv (reverse complementary): only applied to DNA sequence. 
		
		Each row contains a single sequence and each sequence has the same length.
		Lowercase in sequence is automatically converted into uppercase.

		Input example (test.sites):
			GAGGTAAAC
			TCCGTAAGT
			CAGGTTGGA
			ACAGTCAGT
			TAGGTCATT
			TAGGTACTG
			ATGGTAACT
			CAGGTATAC
			TGTGTGAGT
			AAGGTAAGT	
		'''
		if dna:
			self.seq_type = 'DNA'
		else:
			self.seq_type = 'PROTEIN'
		if name is None:
			self.motif_name = 'Unknown'
		else:
			self.motif_name = name
		if rv:
			tab = string.maketrans('ACGT','TGCA')
		self.seq_count = 0.0
		self.data = defaultdict(dict)	#base_position (column of .sites): base_type : base_count
		self.raw_data = defaultdict(list)	#base_position: list of ACGT in each column
		self.seq_lengths = set()
		self.DNA_bases = ['A','C','G','T']
		self.protein_bases = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
		
		for l in open(sites,'r'):
			if l.startswith('#'):continue
			if l.startswith('>'):continue
			l = l.strip(' \r\n').upper()
			
			# check if all bases are valid symbols
			skip = False
			if self.seq_type == 'DNA':
				for b in l:
					if b not in  self.DNA_bases:
						print("Uncognize DNA base: \"%s\" in %s. Skipped." % (b,l), file=sys.stderr)
						skip = True
						break
			elif self.seq_type == 'PROTEIN':
				for b in l:
					if b not in  self.protein_bases:
						print("Uncognize DNA base: \"%s\" in %s. Skipped." % (b,l), file=sys.stderr)
						skip = True
						break
			if skip:
				continue			
			
			self.seq_lengths.add(len(l))
			self.seq_count += 1
			
			if rv:
				l = l.translate(tab)[::-1]
			
			# read sites into dict(dict)
			for i,v in enumerate(l):
				self.raw_data[i].append(v)
				if v not in self.data[i]:
					self.data[i][v] = 1.0
				else:
					self.data[i][v] += 1.0
		
		# check if all sequences have the same length
		if len(self.seq_lengths) != 1:
			print("Sequence lengths are not equal!", file=sys.stderr)
			sys.exit(1)
		else:
			self.motif_length =self.seq_lengths.pop()
	
	def motif_length(self):
		'''
		Return the motif length (nt)
		'''
		return self.motif_length	
	
	def toPFM(self,FOUT=sys.stdout):
		'''
		Convert motif sites data into position frequency matrix (PFM)	
		'''	
		pfm = []
		if self.seq_type == 'DNA':
			bases = self.DNA_bases
		else:
			bases = self.protein_bases		
			
		for i in range(self.motif_length):	# i is motif position starting from 0
			tmp = []	# list of base_count for each column of sites file
			for b in bases:
				if b in self.data[i]:
					tmp.append(self.data[i][b])
				else:
					tmp.append(0.0)
			pfm.append(tmp)
			
		pfm = np.transpose(np.array(pfm))

		#print("\n\n# PSSM matrix", file=FOUT)
		print('Base\t' + '\t'.join([str(i+1) for i in range(self.motif_length)]), file=FOUT)
		for i,name in enumerate(bases):
			print(name + '\t' + '\t'.join([str(j) for j in pfm[i]]), file=FOUT)

	def toJaspar(self,FOUT=sys.stdout):
		'''
		Convert motif sites data into Jaspar format (.pfm)
		
		Jaspar format example:
			> Mycn
			A [ 0 29 0 2 0 0 ]
			C [31 0 30 1 3 0 ]
			G [ 0 0 0 28 0 31]
			T [ 0 2 1 0 28 0 ]		
		'''	
		pfm = []
		if self.seq_type == 'DNA':
			bases = self.DNA_bases
		else:
			bases = self.protein_bases		
			
		for i in range(self.motif_length):	# i is motif position starting from 0
			tmp = []	# list of base_count for each column of sites file
			for b in bases:
				if b in self.data[i]:
					tmp.append(self.data[i][b])
				else:
					tmp.append(0.0)
			pfm.append(tmp)
			
		pfm = np.transpose(np.array(pfm))
		print('> %s' % self.motif_name, file=FOUT)
		for i,b in enumerate(bases):
			print(b + ' [ ' + ' '.join([str(j) for j in pfm[i]]) + ']', file=FOUT)
	
	def toRawPSSM(self, FOUT=sys.stdout):
		'''
		Convert motif sites data into raw PSSM format (.pfm)
		
		raw PSSM format example:
			>Mync
			0 31 0 0
			29 0 0 2
			0 30 0 1
			2 1 28 0
			0 3 0 28
			0 0 31 0		
		'''	
		pfm = []
		if self.seq_type == 'DNA':
			bases = self.DNA_bases
		else:
			bases = self.protein_bases		
			
		for i in range(self.motif_length):	# i is motif position starting from 0
			tmp = []	# list of base_count for each column of sites file
			for b in bases:
				if b in self.data[i]:
					tmp.append(self.data[i][b])
				else:
					tmp.append(0.0)
			pfm.append(tmp)
			
		pfm = np.array(pfm)
		print('>%s' % self.motif_name, file=FOUT)
		for i in range(self.motif_length):
			print(' '.join([str(j) for j in pfm[i]]), file=FOUT)

	def toMEME(self,pseudocount=0.8, FOUT=sys.stdout):
		'''
		Convert motif sites data into meme's position-specific probability matrix
		
		MEME format example:
			------------------------
			Motif 2 position-specific probability matrix
			------------------------
			letter-probability matrix: alength= 4 w= 6 nsites= 31
			0 31 0 0
			29 0 0 2
			0 30 0 1
			2 1 28 0
			0 3 0 28
			0 0 31 0		
		'''	
		pfm = []
		ppm = []

		if self.seq_type == 'DNA':
			bases = self.DNA_bases
		else:
			bases = self.protein_bases
			
		for i in range(self.motif_length):	# i is motif position starting from 0
			tmp = []	# list of base_count for each column of sites file			
			for b in bases:
				if b in self.data[i]:
					tmp.append(self.data[i][b])
				else:
					tmp.append(0.0)
			pfm.append(tmp)
			
		pfm = np.transpose(np.array(pfm))
		pfm = pfm + pseudocount/4.0
		ppm = pfm/pfm.sum(axis=0)			
		
		ppm = np.transpose(ppm)
		
		print('-'*40, file=FOUT)
		print(self.motif_name + ' position-specific probability matrix', file=FOUT)
		print('-'*40, file=FOUT)
		print('letter-probability matrix: alength= %d w= %d nsites= %d' % (len(bases), self.motif_length, self.seq_count), file=FOUT)
		for i in ppm:
			print(' ' + ' '.join([str(j) for j in i]), file=FOUT)
		
	def toPPM(self,pseudocount=0.8, FOUT=sys.stdout):
		'''
		Convert motif sites data into position probability matrix (PPM)
		Default pseudocount of 0.8 is determined from this paper:
		http://nar.oxfordjournals.org/content/37/3/939.full
		'''	
		pfm = []
		ppm = []

		if self.seq_type == 'DNA':
			bases = self.DNA_bases
		else:
			bases = self.protein_bases
			
		for i in range(self.motif_length):	# i is motif position starting from 0
			tmp = []	# list of base_count for each column of sites file			
			for b in bases:
				if b in self.data[i]:
					tmp.append(self.data[i][b])
				else:
					tmp.append(0.0)
			pfm.append(tmp)
			
		pfm = np.transpose(np.array(pfm))
		pfm = pfm + pseudocount/4.0
		ppm = pfm/pfm.sum(axis=0)			

		#print("\n\n# PPM matrix", file=FOUT)
		print('Base\t' + '\t'.join([str(i+1) for i in range(self.motif_length)]), file=FOUT)
		for i,name in enumerate(bases):
			print(name + '\t' + '\t'.join([str(j) for j in ppm[i]]), file=FOUT)

	def toPWM(self,pseudocount=0.8, bg=None, FOUT=sys.stdout):
		'''
		Convert motif sites data into position weight matrix (PWM)
		PWM is a matrix of log likelihood between sites and background.
		
		Default pseudocount of 0.8 is determined from this paper:
		http://nar.oxfordjournals.org/content/37/3/939.full
		
		if bg is "None", universal background will be used:
		 * DNA:
		 	A = C = G = T = 0.25 (i.e. 1/4)
		 * Protein
		 	A = R = N = ... = V = 0.05 (i.e. 1/20)
		 
		 Otherwise, bg is a dictionary with base as key and the corresponding
		 base frequency as value. eg
		 
		 bg = {'A':0.23, 'C':0.26,'G':0.29,'T':0.22}
		
		'''	
		pwm = []
		pfm = []
		ppm = []
		background = {}
		
		if self.seq_type == 'DNA':
			bases = self.DNA_bases
		else:
			bases = self.protein_bases
		
		#determine background frequency
		if bg is None:
			for b in bases:
				background[b] = 1.0/len(bases)
						
		for i in range(self.motif_length):	# i is motif position starting from 0
			tmp = []	# list of base_count for each column of sites file
			for b in bases:
				if b in self.data[i]:
					tmp.append(self.data[i][b])
				else:
					tmp.append(0.0)
			pfm.append(tmp)
			
		pfm = np.transpose(np.array(pfm))
		pfm = pfm + pseudocount/4.0
		ppm = pfm/pfm.sum(axis=0)		
		
		for i,name in enumerate(bases):
			tmp = [np.log(j/background[name]) for j in ppm[i]]
			pwm.append(tmp)
		
		#print("\n\n# PWM matrix", file=FOUT)
		print('Base\t' + '\t'.join([str(i+1) for i in range(self.motif_length)]), file=FOUT)
		for i,name in enumerate(bases):
			print(name + '\t' + '\t'.join([str(np.log(j/background[name])) for j in ppm[i]]), file=FOUT)


	
		
