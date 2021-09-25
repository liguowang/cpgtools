#!/usr/bin/env python3

"""
Description
-----------
This program performs PCA (principal component analysis) for samples.

Example of input data file
---------------------------
CpG_ID	Sample_01	Sample_02	Sample_03	Sample_04
cg_001	0.831035	0.878022	0.794427	0.880911
cg_002	0.249544	0.209949	0.234294	0.236680
cg_003	0.845065	0.843957	0.840184	0.824286

Example of the input group file
---------------------------
Sample,Group
Sample_01,normal
Sample_02,normal
Sample_03,tumor
Sample_04,tumor

Notes
-----
* Rows with missing values will be removed
* Beta values will be standardized into z scores
* Only the first two components will be visualized
"""


import sys
import subprocess
from optparse import OptionParser
from cpgmodule.utils import *
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.10.0"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def pick_colors(n):
	my_colors = ['#e6194B', '#3cb44b', '#4363d8', '#f58231', '#911eb4', '#42d4f4', '#f032e6', '#bfef45', '#fabebe', '#469990', '#e6beff', '#9A6324', '#fffac8', '#800000', '#aaffc3', '#808000', '#ffd8b1', '#000075', '#a9a9a9','#ffe119']
	if n > len(my_colors):
		print ("Only support 21 different colors", file = sys.stderr)
		sys.exit()
	return my_colors[0:n]

	
def main():
	
	usage="%prog [options]" + "\n"
	parser = OptionParser(usage,version="%prog " + __version__)
	parser.add_option("-i","--input_file",action="store",type="string",dest="input_file",help="Tab-separated data frame file containing beta values with the 1st row containing sample IDs and the 1st column containing CpG IDs.")
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Comma-separated group file defining the biological groups of each sample. Different groups will be colored differently in the PCA plot. Supports a maximum of 20 groups.")
	parser.add_option("-n","--ncomponent",action="store",type='int', dest="n_components", default=2, help="Number of components. default=%default" )
	parser.add_option("-l","--label",action="store_true",default=False,dest="text_label",help="If True, sample ids will be added underneath the data point. default=%default")
	parser.add_option("-c","--char",action="store",type='int', default=1, dest="plot_char",help="Ploting character: 1 = 'dot', 2 = 'circle'. default=%default")
	parser.add_option("-a","--alpha",action="store",type='float', default=0.5, dest="plot_alpha",help="Opacity of dots. default=%default")
	parser.add_option("-x","--loc",action="store",type='int', default=1, dest="legend_location",help="Location of legend panel: 1 = 'topright', 2 = 'bottomright', 3 = 'bottomleft', 4 = 'topleft'. default=%default")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="The prefix of the output file.")
	(options,args)=parser.parse_args()
	
	print ()
	if not (options.input_file):
		print (__doc__)
		parser.print_help()
		sys.exit(101)
	
	if not (options.out_file):
		print (__doc__)
		parser.print_help()
		sys.exit(103)	
	if options.n_components < 2:
		options.n_components = 2
	
	pch = {1:20, 2:1}
	legend_pos = {1:'topright', 2: 'bottomright', 3:'bottomleft', 4:'topleft'}
	printlog("Reading input file: \"%s\" ..." % (options.input_file))
	df1 = pd.read_csv(options.input_file, index_col = 0, sep="\t")
	
	#remove NA and transpose
	df2 = df1.dropna(axis=0, how='any')
	printlog("%d rows with missing values were removed." % (len(df1) - len(df2)))
	#print (df2.head())
	
	printlog("Transposing data frame ...")
	df2 = df2.T
	#print (df2.head()) 
	
	printlog("Standarizing values ...")
	x = df2.values
	x = StandardScaler().fit_transform(x)
	
	printlog("Reading group file: \"%s\" ..." % (options.group_file))
	group = pd.read_csv(options.group_file, index_col=0, header=0,names=['Sample_ID', 'Group_ID'])
	group.index = group.index.map(str)
	
	#check if sample IDs are unique
	if len(group.index) != len(group.index.unique()):
		print ("Sample IDs are not unique", file = sys.stderr)
		sys.exit()	
	
	group_names = group['Group_ID'].unique().tolist()	# a list of unique group names
	color_names = pick_colors(len(group_names))	# a list of unique colors
	group_to_col = dict(zip(group_names, color_names))
	color_list = [group_to_col[g] for g in group['Group_ID']]
	group['Colors'] = color_list
		
	pca = PCA(n_components = options.n_components, random_state = 0)
	principalComponents = pca.fit_transform(x)	
	pca_names = [str(i)+str(j) for i,j in zip(['PC']*options.n_components,range(1,options.n_components+1))]
	principalDf = pd.DataFrame(data = principalComponents, columns = pca_names, index = df2.index)	
	
	finalDf = pd.concat([principalDf, group], axis = 1, sort=False)
	finalDf.index.name = 'Sample_ID'
	
	printlog("Writing PCA results to file: \"%s\" ..." % (options.out_file + '.PCA.tsv'))
	finalDf.to_csv(options.out_file + '.PCA.tsv', sep="\t")
	
	pca_vars = pca.explained_variance_ratio_
	for n,v in zip(pca_names, pca_vars):
		print ("Variance explained by %s : %.4f%%" % (n, v*100))
		
	
	ROUT = open(options.out_file + '.PCA.r','w')
	
	print ('pdf(file=\"%s\", width=8, height=8)' % (options.out_file + '.PCA.pdf'),file=ROUT)
	print ('')
	print ('d = read.table(file=\"%s\", sep="\\t", header=TRUE,  comment.char = "", stringsAsFactors=FALSE)' % (options.out_file + '.PCA.tsv'), file=ROUT)
	print ('attach(d)', file=ROUT)
	if options.plot_alpha:
		print ('library(scales)', file=ROUT)
		print ('plot(PC1, PC2, col = alpha(Colors, %f), pch=%d, cex=1.5, main="PCA 2D map")' % (options.plot_alpha, pch[options.plot_char]), file=ROUT)
	else:
		print ('plot(PC1, PC2, col = Colors, pch=%d, cex=1.2, main="PCA 2D map")' % pch[options.plot_char], file=ROUT)	
	
	if options.text_label:
		print ('text(PC1, PC2, labels=Sample_ID, col = Colors, cex=0.5, pos=1)', file=ROUT)
	
	print ('legend("%s", legend=c(%s), col=c(%s), pch=%d,cex=1)' %  (legend_pos[options.legend_location], ','.join(['"' + str(i) + '"' for i in group_names]), ','.join(['"' + str(group_to_col[i]) + '"' for i in group_names]), pch[options.plot_char]), file=ROUT)
	
	
	print ('dev.off()', file=ROUT)
	ROUT.close()
	
	try:
		subprocess.call("Rscript " + options.out_file + '.PCA.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.PCA.r', file=sys.stderr)
	pass

	
if __name__=='__main__':
	main()	
