#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 13:57:36 2023
"""
import sys
import os
import pandas as pd
# import numpy as np
import logging
from sklearn.impute import KNNImputer


def update_df_row(df_in, fill_value='mean'):
    """

    Fill missing values by **row**

    Parameters
    ----------
    df_in : str
        Input dataframe.
    fill_value : str, optional
        How to fill the missing value.
        Must be one of ["mean", "median", "min", "max"]
        The default is 'mean'.

    Returns
    -------
    Pandas DataFrame.

    """
    for cgID, betas in df_in.iterrows():
        if betas.isnull().values.any() and betas.isnull().sum() < len(betas):
            if fill_value == 'mean':
                df_in.loc[cgID, :] = df_in.fillna(betas.mean())
            elif fill_value == 'median':
                df_in.loc[cgID, :] = df_in.fillna(betas.median())
            elif fill_value == 'min':
                df_in.loc[cgID, :] = df_in.fillna(betas.min())
            elif fill_value == 'max':
                df_in.loc[cgID, :] = df_in.fillna(betas.max())
            else:
                logging.error(
                    'Must be one of ["mean", "median", "min", "max"]')
        else:
            continue
    return df_in


def impute_beta(input_df, method, ref=None, k=None, w='uniform'):
    """
    Parameters
    ----------
    input_df : str
        The input tabular structure file containing DNA methylation data.
        #example of CSV file
        ID_REF,A,B,C,
        cg26928153,0.86007,0.79695,np.nan
        cg16269199,np.nan,0.64148,0.65569
        cg13869341,0.76405,np.nan,0.7059
        ...
        #import pandas as pd
        #import numpy as np
        #dat = {'cg26928153':[0.86007,0.79695,np.nan],
                'cg16269199':[np.nan,0.64148,0.65569],
                'cg13869341':[0.76405,np.nan,0.7059]}
        #df1 = pd.DataFrame(dat)

    method : int, optional
        How to replace missing values. Take one of below:
            -1: Remove CpGs with any missing values.
            0: Fill all missign values with "0.0".
            1: Fill all missign values with "1.0".
            2: Fill the missing values with **column mean**
            3: Fill the missing values with **column median**
            4: Fill the missing values with **column min**
            5: Fill the missing values with **column max**
            6: Fill the missing values with **row mean** 
            7: Fill the missing values with **row median**
            8: Fill the missing values with **row min**
            9: Fill the missing values with **row max**
            10: Fill the missing values with **external reference**
            11: Fill the missing values with average values from K Nearest Neighbors (KNN). (default)
    ref : str
        Tab or comma separated file. The first column is CpG ID, the 2nd
        column is beta value.
    k : int
        Number of neighboring samples to use for imputation. If k is None,
        k = sqrt(n) where n is the "number of samples" in the input_df.
        Only effective when used imputation method is "KNN".
    w : str
        Weight function used in prediction. 
        (https://scikit-learn.org/stable/modules/generated/sklearn.impute.KNNImputer.html)
        Only effective when used imputation method is "KNN".

    Returns
    -------
    DataFrame with missing values filled.
    """
    if ref is not None:
        if method != 10:
            logging.warning(
                "Ignore the external reference file \"%s\". \
                To use external file, you must set 'method = 10'." % ref)

    if method == -1:
        logging.info("Remove CpGs with any missing values ...")
        output_df = input_df.dropna(axis=0, how='any')
    elif method == 0:
        logging.info("Fill missing values with ZEROs ...")
        output_df = input_df.fillna(0)
    elif method == 1:
        logging.info("Fill missing values with ONEs ...")
        output_df = input_df.fillna(1)
    elif method == 2:
        logging.info("Fill missing values with column (sample) MEAN ...")
        output_df = input_df.fillna(input_df.mean())
    elif method == 3:
        logging.info("Fill missing values with column (sample) MEDIAN ...")
        output_df = input_df.fillna(input_df.median())
    elif method == 4:
        logging.info("Fill missing values with column (sample) MIN ...")
        output_df = input_df.fillna(input_df.min())
    elif method == 5:
        logging.info("Fill missing values with column (sample) MAX ...")
        output_df = input_df.fillna(input_df.max())
    elif method == 6:
        logging.info("Fill missing values with row (probe) MEAN ...")
        output_df = update_df_row(input_df, 'mean')
    elif method == 7:
        logging.info("Fill missing values with row (probe) MEDIAN ...")
        output_df = update_df_row(input_df, 'median')
    elif method == 8:
        logging.info("Fill missing values with row (probe) MIN ...")
        output_df = update_df_row(input_df, 'min')
    elif method == 9:
        logging.info("Fill missing values with row (probe) MAX ...")
        output_df = update_df_row(input_df, 'max')
    elif method == 10:
        logging.info("Fill missing values with external reference ...")
        if ref is not None and os.path.exists(ref):
            # dict: key is cgID, value is beta
            external_ref = {}
            # read external file
            with open(ref) as external:
                line = external.read().strip()
                if ref.endswith('.csv'):
                    f = line.split(',')
                else:
                    f = line.split()
                try:
                    external_ref[f[0]] = float(f[1])
                except ValueError:
                    pass
            for cgID, betas in input_df.iterrows():
                if betas.isnull().values.any():
                    if cgID in external_ref:
                        ref_beta = external_ref[cgID]
                        input_df.loc[cgID, :] = betas.fillna(ref_beta)
                    else:
                        continue
                else:
                    continue
            output_df = input_df
        else:
            logging.error("External file %s does not exist." % ref)
            sys.exit()
    elif method == 11:
        if k is None:
            nb = int(input_df.shape[1]**0.5)
        else:
            nb = k
        imputer = KNNImputer(n_neighbors=nb, weights=w)
        after = imputer.fit_transform(input_df)
        output_df = pd.DataFrame(after, index = input_df.index, columns = input_df.columns)
    return output_df


if __name__ == '__main__':
    input_df = pd.read_csv(sys.argv[1], index_col=0)
    print(input_df)
    print()
    output_df = impute_beta(input_df, method=10, ref='ref')
    print(output_df)