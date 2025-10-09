#!/usr/bin/env python3

"""
Description
-----------
This program corrects batch effect using the "combat" algorithm:

W. Evan Johnson, et al, Adjusting batch effects in microarray expression data using empirical Bayes methods, Biostatistics, 2007.

Example of input data file
---------------------------
CpG_ID    Sample_01    Sample_02    Sample_03    Sample_04
cg_001    0.831035    0.878022    0.794427    0.880911
cg_002    0.249544    0.209949    0.234294    0.236680
cg_003    0.845065    0.843957    0.840184    0.824286
...

Example of batch file
-------------------------------
Sample,Group
Sample_01,plate_1
Sample_02,plate_1
Sample_03,plate_2
Sample_04,plate_2
...

Notes
-----
* Rows with missing values will be removed

"""


import sys
import subprocess
from optparse import OptionParser
from cpgmodule.utils import *
from cpgmodule._version import __version__
import pandas as pd
#from sklearn.preprocessing import StandardScaler
#from sklearn.decomposition import PCA
from sklearn.impute import KNNImputer
from collections import Counter
from combat.pycombat import pycombat
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from impyutelib import nan_indices


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

def pick_colors(n):
    my_colors = list(mcolors.CSS4_COLORS.keys())
    if n > len(my_colors):
        print ("Too many colors requested", file = sys.stderr)
        sys.exit()
    return my_colors[0:n]

def box_plot(df, s_colors,  out_png, ylab="Beta values", title=""):
    s_names = df.columns
    fig, ax = plt.subplots()
    bplot = ax.boxplot(df, patch_artist=True, tick_labels = s_names)
    for patch, color in zip(bplot['boxes'], s_colors):
        patch.set_facecolor(color)
    ax.set_xticklabels(s_names, rotation='vertical')
    plt.ylabel(ylab)
    plt.title(title)
    plt.savefig(out_png)

def main():
    
    usage="%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input_file",action="store",type="string",dest="input_file",help="Tab-separated data frame file containing beta values with the 1st row containing sample IDs and the 1st column containing CpG IDs.")
    parser.add_option("-k",action="store",type="int", default=3, dest="n_neighbors",help="Number of neighbors to use for imputation. default=%default")
    parser.add_option("--axis",type="choice",choices=[0,1],default=1,dest="axis_choice",help="How to do imputation (using the KNN algorithm) if the input file has missing values. 1: search columns for k nearest neighbours; 0: Search rows for k nearest neighbours. default=%default")
    parser.add_option("-g","--group",action="store",type="string",dest="group_file",help="Comma-separated group file defining the batch groups of each sample.")
    parser.add_option("-o","--output",action="store",type='string', dest="out_file",help="The prefix of the output file.")
    (options,args)=parser.parse_args()
    
    print ()
    if not (options.input_file):
        print (__doc__)
        parser.print_help()
        sys.exit(101)
    if not (options.group_file):
        print (__doc__)
        parser.print_help()
        sys.exit(101)    
    if not (options.out_file):
        print (__doc__)
        parser.print_help()
        sys.exit(103)    

    beta_out1 = options.out_file + '.combat.tsv'
    beta_out2 = options.out_file + '.combat_withNAs.tsv'
    beta_boxplot_before = options.out_file + '.boxplot.png'
    beta_boxplot_after = options.out_file + '.boxplot_combat.png'
    
    printlog("Reading input file: \"%s\" ..." % (options.input_file))
    df1 = pd.read_csv(options.input_file, index_col = 0, sep="\t")
    
    # count NA
    total_na = df1.isna().sum().sum()
    printlog("The input file: \"%s\" contains %d missing values" % (options.input_file, total_na))

    # if there are any NAs, either impute or remove
    if total_na > 0:
        printlog("Imputing missing values using KNN ...")
        na_locations = nan_indices(df1.to_numpy())
        imputer = KNNImputer(n_neighbors = options.n_neighbors)
        if options.axis_choice == 1:
            input_df = df1.T
            after = imputer.fit_transform(input_df)
            df2 = pd.DataFrame(after, index = input_df.index,
                                  columns = input_df.columns).T
            #output_df = output_df.round(args.decimal)
        elif options.axis_choice == 0:
            input_df = df1
            after = imputer.fit_transform(input_df)
            df2 = pd.DataFrame(after, index = input_df.index,
                                  columns = input_df.columns)
        else:
            raise ValueError("axis only accepts 0 or 1.")
    else:
        df2 = df1

    printlog("Reading group file: \"%s\" ..." % (options.group_file))
    group = pd.read_csv(options.group_file, index_col=0, names=['Sample_ID', 'Group_ID'])
    #check if sample IDs are unique
    if len(group.index) != len(group.index.unique()):
        print ("Sample IDs are not unique", file = sys.stderr)
        sys.exit()    
    group.index = group.index.map(str)
    printlog("Group/batch \"%s\" contains %d samples" % (options.group_file, len(group.index)))
    
    # a list of unique group names, and their frequencies
    group_info = Counter(group['Group_ID'])
    print(list(group['Group_ID']))
    # a list of unique colors
    color_names = pick_colors(len(group_info))    
    color_list = []
    for name,count in zip(color_names, list(group_info.values())):
        color_list.extend([name]*count)

    printlog("Generate boxplot before correction. Save to '%s'" % beta_boxplot_before)
    box_plot(df2, s_colors=color_list, out_png=beta_boxplot_before, title="Before batch effects correction")

    # remove batch effect
    printlog("Removing batch effect ...")
    df_corrected = pycombat(df2, list(group['Group_ID']))
    
    printlog("Generate boxplot after correction. Save to '%s'" % beta_boxplot_after)
    box_plot(df_corrected, s_colors=color_list, out_png=beta_boxplot_after, title="After batch effects correction")

    printlog("Save data after combat to \"%s\", keep imputed values ..." % beta_out1)
    df_corrected.to_csv(beta_out1,sep="\t")
    
    # add the original NAs (if any) back.
    # Only perform combat, no KNN imputation
    if total_na > 0:
        for i,j in na_locations:
            df_corrected.iat[i, j] = np.nan
        printlog("Save data after combat to \"%s\", keep orignal missing values ..." % beta_out2)
        df_corrected.to_csv(beta_out2,sep="\t",na_rep='NA')

if __name__=='__main__':
    main()    
