import numpy as np
from impyute.ops import matrix
from impyute.ops import wrapper

@wrapper.wrappers
@wrapper.checks
def mean(data):
    """ Substitute missing values with the mean of that column.

    Parameters
    ----------
    data: numpy.ndarray
        Data to impute.

    Returns
    -------
    numpy.ndarray
        Imputed data.

    """
    nan_xy = matrix.nan_indices(data)
    for x_i, y_i in nan_xy:
        row_wo_nan = data[:, [y_i]][~np.isnan(data[:, [y_i]])]
        new_value = np.mean(row_wo_nan)
        data[x_i][y_i] = new_value
    return data

@wrapper.wrappers
@wrapper.checks
def median(data):
    """ Substitute missing values with the median of that column(middle).

    Parameters
    ----------
    data: numpy.ndarray
        Data to impute.

    Returns
    -------
    numpy.ndarray
        Imputed data.

    """
    nan_xy = matrix.nan_indices(data)
    cols_missing = set(nan_xy.T[1])
    medians = {}
    for y_i in cols_missing:
        cols_wo_nan = data[:, [y_i]][~np.isnan(data[:, [y_i]])]
        median_y = np.median(cols_wo_nan)
        medians[str(y_i)] = median_y
    for x_i, y_i in nan_xy:
        data[x_i][y_i] = medians[str(y_i)]
    return data

@wrapper.wrappers
@wrapper.checks
def mode(data):
    """ Substitute missing values with the mode of that column(most frequent).

    In the case that there is a tie (there are multiple, most frequent values)
    for a column randomly pick one of them.

    Parameters
    ----------
    data: numpy.ndarray
        Data to impute.

    Returns
    -------
    numpy.ndarray
        Imputed data.

    """
    nan_xy = matrix.nan_indices(data)
    modes = []
    for y_i in range(np.shape(data)[1]):
        unique_counts = np.unique(data[:, [y_i]], return_counts=True)
        max_count = np.max(unique_counts[1])
        mode_y = [unique for unique, count in np.transpose(unique_counts)
                  if count == max_count and not np.isnan(unique)]
        modes.append(mode_y)  # Appends index of column and column modes
    for x_i, y_i in nan_xy:
        data[x_i][y_i] = np.random.choice(modes[y_i])
    return data
