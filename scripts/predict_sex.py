#!/usr/bin/env python3

"""
#==============================================================================
Predict sex based on the semi-methylation (also known as genomic imprinting) 
ratio. This method leverages the fact that, due to X chromosome inactivation,
females have a higher proportion of semi-methylated CpGs on their X chromosomes.
A log2(ratio) greater than 0 indicates a female, while a log2(ratio) less than
0 indicates a male.

Example of input data file
---------------------------
CpG_ID    Sample_01    Sample_02    Sample_03    Sample_04
cg_001    0.831035    0.878022    0.794427    0.880911
cg_002    0.249544    0.209949    0.234294    0.236680
cg_003    0.845065    0.843957    0.840184    0.824286

Example of output file
----------------------
Sample_ID    log2_SM_ratio    Predicted_sex
Sample_01    -2.249628052954919      Male
Sample_02    -2.2671726671830674     Male
Sample_03    1.4530581933290616      Female
Sample_04    1.4808015115356654      Female

...

"""
import sys
import numpy as np
from optparse import OptionParser
from cpgmodule.utils import printlog
from cpgmodule import ireader
import pandas as pd
from cpgmodule._version import __version__

__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

    
def main():
    
    usage="%prog [options]" + "\n"
    parser = OptionParser(usage,version="%prog " + __version__)
    parser.add_option("-i","--input_file",action="store", type="string",dest="input_file", help="Tab-separated data frame file containing beta values with the 1st row containing sample IDs and the 1st column containing CpG IDs.")
    parser.add_option("-x","--xprobe",action="store", type="string",dest="xprobe_file", help="File with CpG IDs mapped to the X chromosome, with one probe listed per row.")
    parser.add_option("-c","--cut",action="store", type='float', dest="cutoff", default=0.0, help="The cutoff of log2(SM ratio) to determine the sex prediction. Log2(SM ratio) greater than this cutoff indicates a female, while a log2(ratio) less than this cutoff indicates a male. default=%default")
    parser.add_option("-o","--output",action="store", type='string', dest="out_file", help="The prefix of the output file.")
    (options,args)=parser.parse_args()
    
    print ()
    if not (options.input_file):
        print (__doc__)
        parser.print_help()
        sys.exit(101)
    if not (options.xprobe_file):
        print (__doc__)
        parser.print_help()
        sys.exit(102)    
    if not (options.out_file):
        print (__doc__)
        parser.print_help()
        sys.exit(103)    
    
    printlog("Reading X probes from: \"%s\"" % (options.xprobe_file))
    x_cpgs = set()
    for l in ireader.reader(options.xprobe_file):
        l = l.strip()
        if l.startswith('#'):
            continue
        x_cpgs.add(l)
    printlog("Total %d X probes loaded." % len(x_cpgs))
    
    printlog("Reading input file: \"%s\"" % (options.input_file))
    df1 = pd.read_csv(options.input_file, index_col = 0, sep="\t")
    #print (df1)
    
    #remove any rows with NAs
    df2 = df1.dropna(axis=0, how='any')
    printlog("%d CpGs with missing values were removed." % (len(df1) - len(df2)))
    #print (df2)
    
    sample_cpg_ids = df2.index
    sample_names = df2.columns
    found_x_cpgs = list(x_cpgs & set(sample_cpg_ids))
    printlog("Found %d CpGs located on the chrX from file: %s" % (len(found_x_cpgs), options.input_file))
    
    # only X probes in df3
    df3 = df2.loc[list(found_x_cpgs)]
    #pd.DataFrame.to_csv(df3, options.out_file + '.tmp.tsv', sep="\t", index_label="sample")
    
    low_beta_range = [0, 0.2]
    mid_beta_range = [0.3, 0.7]
    high_beta_range = [0.8, 1.0]
    
    output = {}
    for s in sample_names:
        output[s] = {}
        low_beta_count = pd.cut(df3[s], low_beta_range).count()
        mid_beta_count = pd.cut(df3[s], mid_beta_range).count()
        high_beta_count = pd.cut(df3[s], high_beta_range).count()
        try:
            ratio = np.log2(mid_beta_count/(low_beta_count + high_beta_count))
        except:
           ratio = np.nan
        output[s]['log2_SM_ratio'] = ratio
        
        if ratio > options.cutoff:
            output[s]['Predicted_sex'] = 'Female'
        elif ratio < options.cutoff:
            output[s]['Predicted_sex'] = 'Male'
        else:
            output[s]['Predicted_sex'] = 'Unknown'
    df_out = pd.DataFrame(output).T
    
    outfile = options.out_file + '.predicted_sex.tsv'
    printlog("Writing to file: \"%s\"" % outfile)
    pd.DataFrame.to_csv(df_out, outfile, sep="\t", index_label="Sample_ID")

if __name__=='__main__':
    main()    