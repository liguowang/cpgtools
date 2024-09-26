""" Random utility functions """
from functools import wraps
import numpy as np
import pandas as pd

# Things that get exposed from * import
__all__ = [
    "constantly", "complement", "identity", "thread",
    "execute_fn_with_args_and_or_kwargs", "toy_df",
    "insert_na",
    ]

def thread(arg, *fns):
    if len(fns) > 0:
        return thread(fns[0](arg), *fns[1:])
    else:
        return arg

def identity(x):
    return x

def constantly(x):
    """ Returns a function that takes any args and returns x """
    def func(*args, **kwargs):
        return x
    return func

def complement(fn):
    """ Return fn that outputs the opposite truth values of the
    input function
    """
    @wraps(fn)
    def wrapper(*args, **kwargs):
        return not fn(*args, **kwargs)
    return wrapper

def execute_fn_with_args_and_or_kwargs(fn, args, kwargs):
    """ If args + kwargs aren't accepted only args are passed in"""
    try:
        return fn(*args, **kwargs)
    except TypeError:
        return fn(*args)

def toy_df(nrow, ncol, n_miss, sample_prefix, seed):
    """
    Make a dataFrame (nrow x ncol) with random values between 0 and 1, add
    some missing values (n_miss). Generate a toy dataframe for testing purposes.
    """
    np.random.seed(seed)
    data = np.random.rand(nrow*ncol).reshape((nrow, ncol)).astype(float)
    x_ind = np.random.choice(nrow, n_miss)
    y_ind = np.random.choice(ncol, n_miss)
    for x,y in zip(x_ind, y_ind):
        data[x][y] =  np.nan
    colNames = [sample_prefix + '_' + str(i) for i in range(0,ncol)]
    df = pd.DataFrame(data, columns=colNames)
    return df

def insert_na(df, n_miss, seed):
    np.random.seed(seed)
    nrow,ncol = df.shape
    na_count = 0
    if n_miss >= nrow*ncol:
        out_df = df.replace(df.values, np.nan)
    else:
        tmp = df.to_numpy()
        while(1):
            if na_count >= n_miss:
                break
            x_ind = np.random.choice(nrow)
            y_ind = np.random.choice(ncol)
            if not np.isnan(tmp[x_ind][y_ind]):
                tmp[x_ind][y_ind] = np.nan
                na_count += 1
        out_df = pd.DataFrame(tmp, index=df.index, columns=df.columns)
    return out_df