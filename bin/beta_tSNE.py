#!/usr/bin/env python3

"""
Description
-----------
This program performs t-SNE (t-Distributed Stochastic Neighbor Embedding) analysis for samples.

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
* Different perplexity values can result in significantly different results
"""


import sys
import subprocess
from optparse import OptionParser
from cpgmodule.utils import *
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE

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
	parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Comma-separated group file defining the biological groups of each sample. Different groups will be colored differently in the t-SNE plot. Supports a maximum of 20 groups.")
	parser.add_option("-p","--perplexity",action="store",type='int', dest="perplexity_value", default=5, help="This is a tunable parameter of t-SNE, and has a profound effect on the resulting 2D map. Consider selecting a value between 5 and 50, and the selected value should be smaller than the number of samples (i.e., number of points on the t-SNE 2D map). Default = %default" )
	parser.add_option("-n","--ncomponent",action="store",type='int', dest="n_components", default=2, help="Number of components. default=%default" )
	parser.add_option("--n_iter",action="store",type='int', dest="n_iterations", default=5000, help="The maximum number of iterations for the optimization. Should be at least 250. default=%default" )
	parser.add_option("--learning_rate",action="store",type='float', dest="learning_rate", default=200.0, help="The learning rate for t-SNE is usually in the range [10.0, 1000.0]. If the learning rate is too high, the data may look like a ‘ball’ with any point approximately equidistant from its nearest neighbors. If the learning rate is too low, most points may look compressed in a dense cloud with few outliers. If the cost function gets stuck in a bad local minimum increasing the learning rate may help. default=%default" )
	parser.add_option("-l","--label",action="store_true",default=False,dest="text_label",help="If True, sample ids will be added underneath the data point. default=%default")
	parser.add_option("-c","--char",action="store",type='int', default=1, dest="plot_char",help="Ploting character: 1 = 'dot', 2 = 'circle'. default=%default")
	parser.add_option("-a","--alpha",action="store",type='float', default=0.5, dest="plot_alpha",help="Opacity of dots. default=%default")
	parser.add_option("-x","--loc",action="store",type='int', default=1, dest="legend_location",help="Location of legend panel: 1 = 'topright', 2 = 'bottomright', 3 = 'bottomleft', 4 = 'topleft'. default=%default")
	parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="The prefix of the output file.")
	
	(options,args)=parser.parse_args()
	
	#print (options.text_label)
	#sys.exit(0)
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
	
	n_samples = df1.shape[1]
	#print (n_samples)
	if (options.perplexity_value > n_samples):
		options.perplexity_value = int(n_samples/2)
		printlog("Perplexigty value is set to %d" % options.perplexity_value)
	
	#remove NA and transpose
	df2 = df1.dropna(axis=0, how='any')
	printlog("%d rows with missing values were removed." % (len(df1) - len(df2)))
	#print (df2.head())
	
	printlog("Transposing data frame ...")
	df2 = df2.T
	#print (df2.index) 
		
	printlog("Standarizing values ...")
	x = df2.values
	x = StandardScaler().fit_transform(x)
	#print (x.shape)
	
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
	
	
	tsne = TSNE(n_components = options.n_components, random_state = 0, perplexity = options.perplexity_value, learning_rate = options.learning_rate, n_iter = options.n_iterations)
	tsne_components = tsne.fit_transform(x)
	pc_names = [str(i)+str(j) for i,j in zip(['PC']*options.n_components,range(1,options.n_components+1))]
	principalDf = pd.DataFrame(data = tsne_components, columns = pc_names, index = df2.index)
	principalDf.index.name = 'Sample_ID'
	
	finalDf = pd.concat([principalDf, group], axis=1,sort=False)
	finalDf.index.name = 'Sample_ID'
	
	printlog("Writing t-SNE results to file: \"%s\" ..." % (options.out_file + '.t-SNE.tsv'))
	finalDf.to_csv(options.out_file + '.t-SNE.tsv', sep="\t")
	
		
	ROUT = open(options.out_file + '.t-SNE.r','w')
	
	print ('pdf(file=\"%s\", width=8, height=8)' % (options.out_file + '.t-SNE.pdf'),file=ROUT)
	print ('')
	print ('d = read.table(file=\"%s\", sep="\\t", header=TRUE, comment.char = "", stringsAsFactors=FALSE)' % (options.out_file + '.t-SNE.tsv'), file=ROUT)
	print ('attach(d)', file=ROUT)
	
	if options.plot_alpha:
		print ('library(scales)', file=ROUT)
		print ('plot(PC1, PC2, col = alpha(Colors, %f), pch=%d, cex=1.5, main="t-SNE 2D map")' % (options.plot_alpha, pch[options.plot_char]), file=ROUT)
	else:
		print ('plot(PC1, PC2, col = Colors, pch=%d, cex=1.2, main="t-SNE 2D map")' % pch[options.plot_char], file=ROUT)	
		#print ('plot(PC1, PC2, col = Colors, pch=%d, cex=1, main="t-SNE 2D map")' % pch[options.plot_char], file=ROUT)
	if options.text_label:
		print ('text(PC1, PC2, labels=Sample_ID, col = Colors, cex=0.5, pos=1)', file=ROUT)
	print ('legend("%s", legend=c(%s), col=c(%s), pch=%d,cex=1)' %  (legend_pos[options.legend_location], ','.join(['"' + str(i) + '"' for i in group_names]), ','.join(['"' + str(group_to_col[i]) + '"' for i in group_names]), pch[options.plot_char]), file=ROUT)
	
	
	print ('dev.off()', file=ROUT)
	ROUT.close()
	
	try:
		subprocess.call("Rscript " + options.out_file + '.t-SNE.r', shell=True)
	except:
		print ("Cannot generate pdf file from " + options.out_file + '.t-SNE.r', file=sys.stderr)
	pass

	
if __name__=='__main__':
	main()	
