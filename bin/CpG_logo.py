#!/usr/bin/env python
"""
Description
-----------
This program generates DNA motif logo for a given set of CpGs.
"""


import sys,os
import collections
import subprocess
import numpy as np
import pysam
from optparse import OptionParser
from cpgmodule import ireader
from cpgmodule.utils import *
from cpgmodule import BED
from cpgmodule.imotif import PSSM

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def main():
	print (__doc__)
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input-file",action="store",type="string",dest="input_file",help="BED file specifying the C position. This BED file should have at least 6 columns (Chrom, ChromStart, ChromeEnd, name, score, strand).  Note: Must provide correct *strand* information. This file can be a regular text file or compressed file (.gz, .bz2).")
	parser.add_option("-r","--refgenome",action="store",type="string",dest="genome_file",help="Reference genome seqeunces in FASTA format. Must be indexed using samtools \"faidx\" command. ")
	parser.add_option("-e","--extend",action="store",type="int",dest="extend_size",default=5,help="Number of bases extended to up- and down-stream. default=%default (bp)")
	parser.add_option("-n","--name",action="store",type='string', dest="motif_name",default='motif', help="Motif name. default=%default")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="Prefix of output file.")
	(options,args)=parser.parse_args()
	
	print ()

	if not (options.input_file):
		parser.print_help()
		sys.exit(101)

	if not (options.genome_file):
		parser.print_help()
		sys.exit(102)
	#index refegenome file if it hasn't been done
	if not os.path.exists(options.genome_file + '.fai'):
		printlog("Creating index for %s" % options.genome_file)
		pysam.faidx(options.genome_file)				
	
	if not (options.out_file):
		parser.print_help()
		sys.exit(103)	
	
	refFasta = pysam.Fastafile(options.genome_file)
	FOUT = open(options.out_file + '.fa','w')
	
	printlog("Reading %s ..." % options.input_file)
	for l in ireader.reader(options.input_file):
		if l.startswith('#'):
			continue
		if l.startswith('track'):
			continue
		if l.startswith('browser'):
			continue
		f = l.split()
		if '-' in f:
			strand = '-'
		else:
			strand = '+'
		try:
			chrom = f[0]
			position = int(f[2])
		except:
			print ("BED has at lesat 4 columns. Skip: " + l, file=sys.stderr)	
		
		start = position -  options.extend_size - 1
		end = position + options.extend_size
		if start < 0 or start > end:
			continue
		
		fa_name = '>' + '_'.join([str(i) for i in (chrom,start,end,strand)])
		fa_seq = refFasta.fetch(chrom, start, end).upper()  
		if strand == '-':
			fa_seq = revcomp(fa_seq)
		print (fa_name,file=FOUT)
		print (fa_seq,file=FOUT)
	FOUT.close()
	
	printlog("Generate motif logo ... ")
	try:
		subprocess.call("weblogo  --format PDF -D fasta -c classic -s large -f %s -o %s -t %s" % (options.out_file + '.fa', options.out_file + '.logo.pdf', options.motif_name), shell=True)
		subprocess.call("weblogo  --format PNG -D fasta -c classic -s large -f %s -o %s -t %s" % (options.out_file + '.fa', options.out_file + '.logo.png', options.motif_name), shell=True)
	except:
		print ("Cannot run weblogo. Please install weblogo (https://github.com/WebLogo/weblogo)", file=sys.stderr)
		pass 
	printlog("Motif logo saved to \"%s\" and \"%s\"" % (options.out_file + '.logo.pdf',  options.out_file + '.logo.png'))
	
	
	m = PSSM(sites=options.out_file + '.fa', name = options.motif_name)
	
	printlog("Write position frequency matrix (PFM) to \"%s\"" % (options.out_file + '.pfm'))
	FF = open(options.out_file + '.pfm', 'w')
	m.toPFM(FOUT=FF)
	FF.close()

	printlog("Write position probability matrix (PPM) to \"%s\"" % (options.out_file + '.ppm'))
	FF = open(options.out_file + '.ppm', 'w')
	m.toPPM(FOUT=FF)
	FF.close()

	printlog("Write position weight matrix (PWM) to \"%s\"" % (options.out_file + '.pwm'))
	FF = open(options.out_file + '.pwm', 'w')
	m.toPWM(FOUT=FF)
	FF.close()
		
	printlog("Write Jaspar format matrix to \"%s\"" % (options.out_file + '.jaspar'))
	FF = open(options.out_file + '.jaspar', 'w')
	m.toJaspar(FOUT=FF)
	FF.close()	

	printlog("Write MEME format matrix to \"%s\"" % (options.out_file + '.meme'))
	FF = open(options.out_file + '.meme', 'w')
	m.toMEME(FOUT=FF)
	FF.close()	
	
    	    
if __name__=='__main__':
	main()		
		