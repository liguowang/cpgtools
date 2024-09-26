#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 10:55:14 2022

@author: Liguo Wang
"""
from cpgmodule import ireader
#import sys,os

class MethylSig():
	"""
	Pack DNA methylation signature file into object.

	>>> from cpgmodule import methylClock
	>>> a = methylClock.MethylAge(signature_file = 'coefBlup.tsv', signature_name = 'BLUP', signature_info="")
	>>> a.name
	'BLUP'
	>>> a.Intercept
	91.15396
	>>> a.ncpg
	319607
	"""

	def __init__(self, signature_file, signature_name, tissues = [], unit = '', signature_info = '', reference = '', pub_link = '', method = ''):
		self.name = signature_name
		self.info = signature_info
		self.tissues = tissues
		self.unit = unit
		self.coef = {}
		self.cpgs = []
		self.ncpg = 0
		self.Intercept = 0.0
		self.ref = reference
		self.pubmed = pub_link
		self.method = method
		for l in ireader.reader(signature_file):
			if l.startswith('#'):
				continue
			f = l.split()
			if l.startswith('Intercept'):
				try:
					self.Intercept = float(f[1])
				except:
					self.Intercept = 0.0
			else:
				self.cpgs.append(f[0])
				self.ncpg  += 1
				try:
					self.coef[f[0]] = float(f[1])
					#self.ncpg  += 1
				except:
					continue
