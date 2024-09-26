#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 21:14:40 2024

@author: m102324
"""

import sys
import pandas as pd
import numpy as np
import logging
import argparse
from sklearn.impute import KNNImputer
#from sklearn.experimental import enable_iterative_imputer
#from sklearn.impute import IterativeImputer

from impyute.cs.fast_knn import fast_knn
from impyute.cs.em import em
from impyute.ops.util import toy_df,insert_na
from impyute.cs.random import random_impute
from impyute.cs.buck_iterative import buck_iterative
from missingpy import MissForest

#use pip to install fancyimpute
from fancyimpute import NuclearNormMinimization, SoftImpute, BiScaler

from cpgmodule._version import __version__


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"

logging.basicConfig(
    format="%(asctime)s [%(levelname)s]  %(message)s",
    datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
            
def read_df(infile):
    logging.info("Read input file \"%s\" as a data frame ..." % infile)
    df = pd.read_csv(infile, index_col=0, sep=None, engine='python')
    return df

def nafiller():
    """
    Generic and Specific Methods for DNA Methylation Missing Value Imputation.
    """
   
    general_help = "DNA methylation data imputation package"

    # sub commands and help.
    commands = {
        'DropNA': "Drop rows or columns with NAs.",
        'FillValue': "Replace missing values with a given value.",
        'FillMean': "Impute missing values with row-wise or column-wise means..",
        'FillMedian': "Impute missing values with row-wise or column-wise medians.",
        'FillMin': "Impute missing values with row-wise or column-wise minimum values.",
        'FillMax': "Impute missing values with row-wise or column-wise maximum values",
        'FillRand': "Impute missing values with randomly selected values from \
                    the same row or column.",
        'FillRef': "Impute missing values using values from an external \
                    reference dataset.",
        'KNN': "Impute missing values using scikit-learn's KNNImputer function. \
                Note: slow for large datasets.",
        'KNN2': "Impute missing values using KNN2",
        'fKNN': "Impute missing values using Impyute's fast KNN (fKNN) method.",
        'EM': "Impute missing values using the Expectation Maximization (EM) \
                algorithm.",
        'Buck': "Impute missing values using Buck's method.",
        'NNM': "Impute missing values using Nuclear Norm Minimization (NNM).",
        'SoftImpute': "Impute missing values by iterative soft thresholding of SVD decompositions.",
        'RF': "Impute missing values using Random Forest (RF) prediction.",
        'ToyDf': "Generate a toy dataframe with specified missing values for testing.",
        'InsertNA': "Insert n missing values into an exist dataframe.",
         }
    # create parse
    parser = argparse.ArgumentParser(
        description=general_help, epilog='',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument(
        '-v', '--version', action='version', version='%s %s' %
        ('nafiller', __version__)
        )

    # create sub-parser
    sub_parsers = parser.add_subparsers(help='Sub-command description:')
    
    DropNA_parser = sub_parsers.add_parser('DropNA', help=commands['DropNA'])
    FillValue_parser = sub_parsers.add_parser('FillValue', help=commands['FillValue'])
    FillMean_parser = sub_parsers.add_parser('FillMean', help=commands['FillMean'])
    FillMedian_parser = sub_parsers.add_parser('FillMedian', help=commands['FillMedian'])
    FillMin_parser = sub_parsers.add_parser('FillMin', help=commands['FillMin'])
    FillMax_parser = sub_parsers.add_parser('FillMax', help=commands['FillMax'])
    FillRand_parser = sub_parsers.add_parser('FillRand', help=commands['FillRand'])
    FillRef_parser = sub_parsers.add_parser('FillRef', help=commands['FillRef'])
    KNN_parser = sub_parsers.add_parser('KNN', help=commands['KNN'])
    fKNN_parser = sub_parsers.add_parser('fKNN', help=commands['fKNN'])
    EM_parser = sub_parsers.add_parser('EM', help=commands['EM'])
    Buck_parser = sub_parsers.add_parser('Buck', help=commands['Buck'])
    NNM_parser = sub_parsers.add_parser('NNM', help=commands['NNM'])
    SoftImpute_parser = sub_parsers.add_parser('SoftImpute', help=commands['SoftImpute'])
    RF_parser = sub_parsers.add_parser('RF', help=commands['RF'])
    ToyDf_parser = sub_parsers.add_parser('ToyDf', help=commands['ToyDf'])
    InsertNA_parser = sub_parsers.add_parser('InsertNA', help=commands['InsertNA'])
    
    
    DropNA_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    DropNA_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    DropNA_parser.add_argument(
        '-a', '--axis', type=int, choices=range(2), default=0,
        help="0 : drop rows with any missing values, 1 : drop columns with \
            missing values. Default: 0")
    DropNA_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    FillValue_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    FillValue_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    FillValue_parser.add_argument(
        '-s', '--score', type=float, default=0.0,
        help="The value uesd to fill all NAs.")
    FillValue_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    FillMean_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    FillMean_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    FillMean_parser.add_argument(
        '-a', '--axis', type=int, choices=range(2), default=1,
        help="0 means column, 1 means row. Default: fill NAs with row means")
    FillMean_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    FillMedian_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    FillMedian_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    FillMedian_parser.add_argument(
        '-a', '--axis', type=int, choices=range(2), default=1,
        help="0 means column, 1 means row. Default: fill NAs with row medians")
    FillMedian_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    FillMin_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    FillMin_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    FillMin_parser.add_argument(
        '-a', '--axis', type=int, choices=range(2), default=1,
        help="0 means column, 1 means row. Default: fill NAs with the minimum value of the rows.")
    FillMin_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")
    
    FillMax_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    FillMax_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    FillMax_parser.add_argument(
        '-a', '--axis', type=int, choices=range(2), default=1,
        help="0 means column, 1 means row. Default: fill NAs with the maximum value of the rows.")
    FillMax_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")
    
    FillRand_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    FillRand_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    FillRand_parser.add_argument(
        '-a', '--axis', type=int, choices=range(2), default=1,
        help="0 means column, 1 means row. Default: fill NAs with values randomly selected from rows.")
    FillRand_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    FillRef_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    FillRef_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    FillRef_parser.add_argument(
        '-r', '--ref', type=str,
        help="File name of the external reference.")
    FillRef_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    KNN_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    KNN_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    KNN_parser.add_argument(
        '-k', '--neighbours', type=int,
        help="Number of neighboring samples to use for imputation. If k is \
            None, k = sqrt(n) where n is the \"total number of samples\".")
    KNN_parser.add_argument(
        '-w', '--weightfunction', type=str,
        choices=['uniform', 'distance'],
        default='distance',
        help="Weight function used in predictionaction.")
    KNN_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    fKNN_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    fKNN_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    fKNN_parser.add_argument(
        '-k', '--neighbours', type=int,
        help="Number of neighboring samples to use for imputation. If k is \
            None, k = sqrt(n) where n is the \"total number of samples\".")
    fKNN_parser.add_argument(
        '--eps', type=float,
        default=0,
        help="Refer to the docs for [`scipy.spatial.KDTree.query`]. Must be\
            non-negative float number.")
    fKNN_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    EM_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    EM_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    EM_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")
    
    Buck_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    Buck_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    Buck_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")
    
    NNM_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    NNM_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    NNM_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    SoftImpute_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    SoftImpute_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    SoftImpute_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    RF_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    RF_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    RF_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")

    ToyDf_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Name of the output data frame.")
    ToyDf_parser.add_argument(
        '-r', '--nrow', type=int, default=10,
        help="Number of rows. default: %(default)s")
    ToyDf_parser.add_argument(
        '-c', '--ncol', type=int, default=10,
        help="Number of columns. default: %(default)s")
    ToyDf_parser.add_argument(
        '--na', type=int, default=5,
        help="Number of missing values ingested into the dataframe. default: %(default)s")
    ToyDf_parser.add_argument(
        '-s', '--seed', type=int, default=123,
        help="Seed used to initialize a pseudorandom number generator. default: %(default)s")
    ToyDf_parser.add_argument(
        '--prefix', type=str, default='s',
        help="Prefix of the column names, a series numbers will be appended to the prefix. default: %(default)s")
    ToyDf_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")
    
    InsertNA_parser.add_argument(
        'input', type=str, metavar='input_df',
        help="Input data frame.")
    InsertNA_parser.add_argument(
        'output', type=str, metavar='out_df',
        help="Output data frame.")
    InsertNA_parser.add_argument(
        '--na', type=int,
        help="Number of missing values ingested into the dataframe.")
    InsertNA_parser.add_argument(
        '-s', '--seed', type=int, default=123,
        help="Seed used to initialize a pseudorandom number generator. default: %(default)s")
    InsertNA_parser.add_argument(
        '--decimal', type=int, default=5,
        help="Number of decimal places to round each column to. default: %(default)s")
    
    axis_name = {0:'columns', 1:'rows'}
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(0)
    elif len(sys.argv) >= 2:
        command = sys.argv[1]
        if command.lower() == 'dropna':
            input_df = read_df(args.input)
            #reset
            axis_name = {0:'rows', 1:'columns'}
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Remove %s with any missing values ..." % axis_name[args.axis])
            output_df = input_df.dropna(axis=args.axis, how='any').round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fillvalue':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values with \"%f\" ..." % args.score)
            output_df = input_df.fillna(args.score).round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                         (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fillmean':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            if args.axis == 0:
                logging.info("Replace missing values with COLUMN means ...")
                output_df = input_df.fillna(input_df.mean())
            elif args.axis == 1:
                logging.info("Replace missing values with ROW means ...")
                input_df = input_df.T
                output_df = input_df.fillna(input_df.mean()).T
                logging.info("Unknown parameter.")
                pass
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fillmedian':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            if args.axis == 0:
                logging.info("Replace missing values with COLUMN medians ...")
                output_df = input_df.fillna(input_df.median())
            elif args.axis == 1:
                logging.info("Replace missing values with ROW medians ...")
                input_df = input_df.T
                output_df = input_df.fillna(input_df.median()).T
            else:
                logging.info("Unknown parameter.")
                pass
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fillmin':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            if args.axis == 0:
                logging.info("Replace missing values with COLUMN minimums ...")
                output_df = input_df.fillna(input_df.min())
            elif args.axis == 1:
                logging.info("Replace missing values with ROW minimums ...")
                input_df = input_df.T
                output_df = input_df.fillna(input_df.min()).T
            else:
                logging.info("Unknown parameter.")
                pass
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fillmax':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            if args.axis == 0:
                logging.info("Replace missing values with COLUMN maximums ...")
                output_df = input_df.fillna(input_df.max())
            elif args.axis == 1:
                logging.info("Replace missing values with ROW maximums ...")
                input_df = input_df.T
                output_df = input_df.fillna(input_df.max()).T
            else:
                logging.info("Unknown parameter.")
                pass
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fillrand':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            if args.axis == 0:
                logging.info("Replace missing values with values randomly selected  from the same column ...")
                output_df = random_impute(input_df)
            elif args.axis == 1:
                logging.info("Replace missing values with values randomly selected  from the same row ...")
                input_df = input_df.T
                output_df = random_impute(input_df).T
            else:
                logging.info("Unknown parameter.")
                pass
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fillref':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))

            input_df.index = input_df.index.map(str)
            external_ref = {}
            logging.info("Read external reference file \"%s\" ..." % args.ref)
            for l in open(args.ref,'r'):
                l = l.strip()
                f = l.split()
                try:
                    external_ref[f[0]] = float(f[1])
                except ValueError:
                    pass
            logging.info("Replace missing values with values from external reference ...")
            for ID, betas in input_df.iterrows():
                if betas.isnull().values.any():
                    if ID in external_ref:
                        ref_beta = external_ref[ID]
                        input_df.loc[ID, :] = betas.fillna(ref_beta)
                    else:
                        continue
                else:
                    continue
            output_df = input_df
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'knn':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values using scikit-learn's KNNImputer function ...")
            if args.neighbours is None:
                nb = int(input_df.shape[1]**0.5)
            else:
                nb = args.neighbours
            imputer = KNNImputer(n_neighbors=nb, weights=args.weightfunction)
            after = imputer.fit_transform(input_df)
            output_df = pd.DataFrame(after, index = input_df.index,
                                     columns = input_df.columns)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'fknn':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values using Impyute's fast KNN (fKNN) method ...")
            if args.neighbours is None:
                nb = int(input_df.shape[1]**0.5)
            else:
                nb = args.neighbours
            output_df = fast_knn(input_df, k=nb, eps=args.eps, p=2, 
                                 distance_upper_bound=np.inf, leafsize=10)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'em':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values using EM algorithem ...")
            output_df = em(input_df)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                         (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'buck':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values using Buck's method ...")
            output_df = buck_iterative(input_df)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'nnm':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values using nuclear norm minimization ...")
            X_filled = NuclearNormMinimization().fit_transform(input_df.to_numpy())
            output_df = pd.DataFrame(X_filled, index=input_df.index, columns=input_df.columns)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'softimpute':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values using SoftImpute ...")
            X_filled = SoftImpute().fit_transform(input_df.to_numpy())
            output_df = pd.DataFrame(X_filled, index=input_df.index, columns=input_df.columns)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))

        elif command.lower() == 'rf':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.input, input_df.isna().sum().sum()))
            logging.info("Replace missing values using Random Forest ...")
            imputer = MissForest()
            X_filled = imputer.fit_transform(input_df.to_numpy())
            output_df = pd.DataFrame(X_filled, index=input_df.index, columns=input_df.columns)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))

        elif command.lower() == 'toydf':
            logging.info("Generate toy dataframe ...")
            output_df = toy_df(nrow = args.nrow, ncol = args.ncol, 
                               n_miss = args.na, sample_prefix=args.prefix,
                               seed=args.seed)
            #print(output_df)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                        (args.output, output_df.isna().sum().sum()))
        elif command.lower() == 'insertna':
            input_df = read_df(args.input)
            logging.info("File \"%s\" contains %d missing values ..." %
                         (args.input, input_df.isna().sum().sum()))
            logging.info("Insert %d NAs into dataframe ..." % args.na)
            output_df = insert_na(df=input_df, n_miss=args.na, seed=args.seed)
            output_df = output_df.round(args.decimal)
            output_df.to_csv(args.output, sep="\t", na_rep="NaN")
            logging.info("File \"%s\" contains %d missing values ..." %
                         (args.output, output_df.isna().sum().sum()))

    else:
        print("Unknown command!")
        parser.print_help(sys.stderr)
        sys.exit(0)

if __name__=='__main__':
    nafiller()
