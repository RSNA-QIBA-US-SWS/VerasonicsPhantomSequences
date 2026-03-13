"""
    
OutlierRemoval.py

This script contains helper functions to remove outlier from data vectors.

Author: Kaden Bock (04/12/2025)
    
"""

import numpy as np

def remove_outliers_based_on_set_thresholds(input_data:np.array, thresh:float|tuple):
    """ Removes outliers based on set thresholds. Removed values are set to np.nan.
    
        *** Note: *** This function is currently only set up to operate over 1-D arrays.
        
        If only one value is passed into the thresh argument, then the threshold is assumed
        to be an upper limit. Otherwise, the first argument is assumed to be the lower lim,
        and the second argument is assumed to be the upper lim.

    Args:
        input_data (np.array): The data array over which invalid values are removed.
        thresh (float | tuple): The threshold or set of thresholds over or under which
            to remove invalid points in the data.

    Raises:
        KeyError: In the event that more than 2 thresholds are passed into the function.
        ValueError: In the event that the upper limit is higher than the lower limit.

    Returns:
        np.array: The input data array, will all invalid data points changed to nans.
    """
    
    if type(thresh) is tuple:
        if len(thresh) > 2:
            raise KeyError("Only 2 thresholds can be passed into this function")
        
        lower_lim = thresh[0]
        upper_lim = thresh[1]
        
        if upper_lim < lower_lim:
            raise ValueError("The upper threshold must be greater than the lower threshold")
        
    else:
        lower_lim = -1*np.inf
        upper_lim = thresh
        
    input_data[(input_data < lower_lim)|(input_data>upper_lim)] = np.nan
    
    return input_data

def remove_outliers_based_on_std(input_data:np.array, thresh:float, axis_num:int=0):
    """ Removes outliers based on falling above or below a certain number of the std
        of the unfiltered data.

    Args:
        input_data (np.array): The data array over which invalid values are removed.
        thresh (float): The threshold over or under which to remove invalid points
            in the data. Multiplied by the data's std to set the final threshold.
        axis_num (int, optional): Axis over which to remove outliers. Defaults to 0.

    Returns:
        np.array: The input data array, will all invalid data points changed to nans.
    """
    
    start_mean = np.nanmean(input_data, axis=axis_num)
    start_std = np.nanstd(input_data, axis=axis_num)     
    
    input_data[(input_data<(start_mean - thresh*start_std))|(input_data>(thresh*start_std + start_mean))] = np.nan
    
    return input_data