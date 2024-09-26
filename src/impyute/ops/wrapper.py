""" Decorator functions to wrap around entry and exit

... to easily apply to a function, functions that check/process inputs
and outputs
"""
from functools import wraps
import numpy as np
import pandas as pd

from . import error
from . import matrix
from . import util as u

## Hacky way to handle python2 not having `ModuleNotFoundError`
# pylint: disable=redefined-builtin, missing-docstring
try:
    raise ModuleNotFoundError
except NameError:
    class ModuleNotFoundError(Exception):
        pass
except ModuleNotFoundError:
    pass
# pylint: enable=redefined-builtin, missing-docstring


def handle_df(fn):
    """ Decorator to handle pandas Dataframe object as input

    If the first arg is a pandas dataframe, convert it to a numpy array
    otherwise don't do anything. Cast back to a pandas Dataframe after
    the imputation function has run
    """
    @wraps(fn)
    def wrapper(*args, **kwargs):
        is_df = False
        ## convert tuple to list so args can be modified
        args = list(args)
        ## Either make a copy or use a pointer to the original
        if kwargs.get('inplace'):
            args[0] = args[0]
        else:
            args[0] = args[0].copy()

        ## If input data is a dataframe then cast the input to an np.array
        ## and set an indicator flag before continuing
        if isinstance(args[0], pd.DataFrame):
            is_df = True
            in_ind = args[0].index
            in_columns = args[0].columns
            args[0] = args[0].to_numpy()

        ## function invokation
        results = u.execute_fn_with_args_and_or_kwargs(fn, args, kwargs)

        ## cast the output back to a DataFrame.
        if is_df:
            results = pd.DataFrame(results, index=in_ind, columns=in_columns)
        return results
    return wrapper

def add_inplace_option(fn):
    """ Decorator for inplace option

    Functions wrapped by this can have an `inplace` kwarg to use either a copy of
    data or reference """
    @wraps(fn)
    def wrapper(*args, **kwargs):
        """ Run input checks"""
        ## convert tuple to list so args can be modified
        args = list(args)
        ## Either make a copy or use a pointer to the original
        if kwargs.get('inplace'):
            args[0] = args[0]
        else:
            args[0] = args[0].copy()

        ## function invokation
        return u.execute_fn_with_args_and_or_kwargs(fn, args, kwargs)
    return wrapper

def conform_output(fn):
    """ Decorator to handle impossible values

    Adds two optional kwargs, `coerce_fn` and `valid_fn`.

    `valid_fn` function stub

        def my_coerce_fn(some_literal) -> boolean

    `coerce_fn` function stub

        def my_coerce_fn(arr, x_i, y_i) -> some_literal

    Valid function is something run on each element of the, this is
    the function that we use to indicate whether the value is valid
    or not

    Coerce function has three arguments, the original matrix and
    the two indices of the invalid value x_i and y_i. This function
    will be run on all invalid values.
    """
    @wraps(fn)
    def wrapper(*args, **kwargs):
        def raise_error(arr, x_i, y_i):
            raise error.BadOutputError("{} does not conform".format(arr[x_i, y_i]))
        ## convert tuple to list so args can be modified
        args = list(args)
        # function that checks if the value is valid
        valid_fn = kwargs.get("valid_fn", u.constantly(True))
        # function that modifies the invalid value to something valid
        coerce_fn = kwargs.get("coerce_fn", raise_error)

        ## function invokation
        results = u.execute_fn_with_args_and_or_kwargs(fn, args, kwargs)

        # check each value to see if it's valid
        bool_arr = matrix.map_nd(u.complement(valid_fn), results)
        # get indices of invalid values
        invalid_indices = np.argwhere(bool_arr)
        # run the coerce fn on each invalid indice
        for x_i, y_i in invalid_indices:
            results[x_i, y_i] = coerce_fn(results, x_i, y_i)

        return results
    return wrapper

def wrappers(fn):
    """ Helper decorator, all wrapper functions applied to modify input (matrix
    with missing values) and output (matrix with imputed values)

    NOTE: `handle_df` has to be last as it needs to be in the outer loop (first
    entry point) since every other function assumes you're getting an np.array
    as input
    """
    return u.thread(
        fn,                 # function that's getting wrapped
        add_inplace_option, # allow choosing reference/copy
        conform_output,     # allow enforcing of some spec on returned outputs
        handle_df,          # if df type, cast to np.array on in and df on out
    )

def _shape_2d(data):
    """ True if array is 2D"""
    return len(np.shape(data)) == 2

def _shape_3d(data):
    """ True if array is 3D"""
    return len(np.shape(data)) == 3

def _is_ndarray(data):
    """ True if the array is an instance of numpy's ndarray"""
    return isinstance(data, np.ndarray)

def _dtype_float(data):
    """ True if the values in the array are floating point"""
    return data.dtype == float

def _nan_exists(data):
    """ True if there is at least one np.nan in the array"""
    nan_xy = matrix.nan_indices(data)
    return len(nan_xy) > 0

def checks(fn):
    """ Throw exception if error runs"""
    @wraps(fn)
    def wrapper(*args, **kwargs):
        data = args[0]
        if len(np.shape(data)) != 2:
            raise error.BadInputError("No support for arrays that aren't 2D yet.")
        elif not _shape_2d(data):
            raise error.BadInputError("Not a 2D array.")
        elif not _is_ndarray(data):
            raise error.BadInputError("Not a np.ndarray.")
        elif not _dtype_float(data):
            raise error.BadInputError("Data is not float.")
        elif not _nan_exists(data):
            raise error.BadInputError("No NaN's in given data")
        return u.execute_fn_with_args_and_or_kwargs(fn, args, kwargs)
    return wrapper
