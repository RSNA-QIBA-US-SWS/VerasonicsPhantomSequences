"""
    
QIBA.py

This script contains helper functions to create space-time plots at specific
axial depths during the processing of ARFI-SWEI data.

Author: Kaden Bock (04/05/2025)
    
"""

import numpy as np
import scipy.signal as sps
import scipy.interpolate as spi

def generate_disp_plane_at_depth(disp_data:np.array, axial_ax_mm:np.array, depth_to_generate_plane_mm:float, 
                                depth_of_field_to_average_over:float=0, ax_axis_num:int=0):
    """ Generates a space-time plot of displacement data over a set axial depth.
    
    The output space-time plot has lateral position in the 0th axis, and time in 1st axis. 

    Args:
        disp_data (np.array): The data array containing the displacement estimates.
        axial_ax_mm (np.array): The corresponding vector containing axial positions for the data.
        depth_to_generate_plane_mm (float): The axial dept at which to generate the space-time plot.
        depth_of_field_to_average_over (float, optional): The axial depth of field over which to
            average displacement data over. Defaults to 0, in which data is only taken from the
            single axial trace corresponding to the closest axial trace to the input depth.
        ax_axis_num (int, optional): The number of the axis that corresponds to the axial axis in the
            displacement data array. Defaults to 0.
            
    Raises:
        ValueError: Raises a Value Error in the event that the axis vector over which axial positional
            information is stored is not the same length as the dimensions of the axial position in the
            array in the data plane.

    Returns:
        np.array: A lateral position by time array containing the space-time plot at the specified
            axial depth.
    """
    
    if axial_ax_mm.shape[0] != disp_data.shape[ax_axis_num]:
        raise ValueError("The length of the axial axis and the axis over which axial data"+
                         " is stored in  the displacement data array are not the same. Please"+
                         " double check and correct your inputs.")
   
    closest_ax_idx = np.argmin(np.abs(axial_ax_mm - depth_to_generate_plane_mm))
    
    dAx_mm = np.mean(np.diff(axial_ax_mm))
    ax_DOF_avg_num_pixels = int((depth_of_field_to_average_over/dAx_mm)/2)
    
    disp_data_slicer = [slice(disp_data.shape[ind]) for ind in range(len(disp_data.shape))]
    
    start = int(closest_ax_idx-ax_DOF_avg_num_pixels)
    stop = int(closest_ax_idx+ax_DOF_avg_num_pixels)
    
    disp_data_slicer[ax_axis_num] = slice(start, stop)
    
    disp_plane = np.mean(disp_data[tuple(disp_data_slicer)], axis=ax_axis_num)
    
    return disp_plane

def split_disp_plane_down_axis(disp_plane:np.array, ax:np.array, val:float, ax_num:int, flip_left:bool=True):
    """ Splits a space-time plot down the position specified. Generally used to split down the central lateral
        position.
    

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        ax (np.array): The vector containing positioning information (spacial or temporal) corresponding to the
            axis over which the space-time plot is going to be split.
        val (float): The value at which the split should occur.
        ax_num (int): The axis number over which the split should occur.
        flip_left (bool, optional): Toggle switch to determine if data from one side of the split plane
            should be flipped, such that the wave propagates the same in both directions. Defaults to True.

    Raises:
        ValueError: Raises a Value Error in the event that the axis vector over which to split the
            plane is not the same length as the dimensions of the array in the data plane.

    Returns:
        np.array: The "left" lateral position by time array containing the space-time plot at the specified
            axial depth, split at the input specification.
        np.array: The lateral positioning array corresponding to the "left" side of the split plot.
        np.array: The "right" lateral position by time array containing the space-time plot at the specified
            axial depth, split at the input specification.
        np.array: The lateral positioning array corresponding to the "right" side of the split plot.
    """
    
    if ax.shape[0] != disp_plane.shape[ax_num]:
        raise ValueError("The length of the unit axis and the axis over which to"+
                         " limit the extent of the disp plane are not the same. Please"+
                         " double check and correct your inputs.")
    
    closest_lat_idx = np.argmin(np.abs(ax - val))
    
    disp_plane_left_slicer = [slice(disp_plane.shape[ind]) for ind in range(len(disp_plane.shape))]
    disp_plane_right_slicer = [slice(disp_plane.shape[ind]) for ind in range(len(disp_plane.shape))]
    
    disp_plane_left_slicer[ax_num] = slice(0, closest_lat_idx)
    disp_plane_right_slicer[ax_num] = slice(closest_lat_idx, None)
    
    left_disp_plane = disp_plane[tuple(disp_plane_left_slicer)]
    right_disp_plane = disp_plane[tuple(disp_plane_right_slicer)]
    
    left_ax = ax[0:closest_lat_idx]
    right_ax = ax[closest_lat_idx:]
    
    if flip_left:
        left_disp_plane = np.flip(left_disp_plane, axis=ax_num)
        left_ax = -1*np.flip(left_ax)
    
    return left_disp_plane, left_ax, right_disp_plane, right_ax

def limit_extent_of_disp_plane(disp_plane:np.array, limit_ax:np.array, limits:float|tuple,
                               limit_axis_num:int):
    """ Limits a displacement plane over a set of input limits corresponding to the positional vector.

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        limit_ax (np.array): The vector containing positioning (spacial or temporal) corresponding to the
            axis over which the space-time plot is going to be limited.
        limits (float | tuple): The limits over which the extent of the plane is limited. If a single value
            is provided, the single input is assumed to be an upper limit, with no lower limit used to limit
            the extent of the plane.
        limit_axis_num (int): The number of the axis in the data plane over which the plane should be limited.

    Raises:
        ValueError: Raises a Value Error in the event that the axis vector over which to limit the
            plane is not the same length as the dimensions of the array in the data plane, or in the event
            that the length of the input limits is larger than 2.

    Returns:
        np.array: A lateral position by time array containing the space-time plot at the specified
            axial depth, limited by the input specifications.
        np.array: The axis over which limits were applied, adjusted to account for the applied limits.
    """
    
    if limit_ax.shape[0] != disp_plane.shape[limit_axis_num]:
        raise ValueError("The length of the unit axis and the axis over which to"+
                         " limit the extent of the disp plane are not the same. Please"+
                         " double check and correct your inputs.")
        
    if type(limits) is tuple:
        if len(limits) > 2:
            raise ValueError("This function can only limit the extent of the plane over an upper and lower limit.")
        lower_lim = limits[0]
        upper_lim = limits[1]
    else:
        lower_lim = -1*np.inf
        upper_lim = limits
    
    start_idx = np.argmin(np.abs(limit_ax - lower_lim))
    end_idx = np.argmin(np.abs(limit_ax - upper_lim))
    
    disp_plane_slicer = [slice(disp_plane.shape[ind]) for ind in range(len(disp_plane.shape))]
    disp_plane_slicer[limit_axis_num] = slice(start_idx, end_idx+1)
    
    return disp_plane[tuple(disp_plane_slicer)], limit_ax[slice(start_idx, end_idx+1)]

def remove_and_fix_reverb_frames(disp_plane:np.array, t_ms:np.array, n_additional_steps_to_remove:int=0, 
                                 num_refs:int=None, lat_axis_num:int=0, t_axis_num:int=1):
    """ Removes and fixes offsets caused by large displacement traces caused by ARFI push reverberations.

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        t_ms (np.array): The vector containing the temporal information corresponding to the space-time plot.
        n_additional_steps_to_remove (int, optional): Number of additional frames to remove after the reverb frame.
            Only necessary if the push was strong enough to cause many additional reverb frames. Defaults to 0.
        num_refs (int, optional): Number of reference frames taken before the push event. Defaults to None.
        lat_axis_num (int, optional): The number of the axis that corresponds to the lateral axis in the
            displacement data array. Defaults to 0.
        t_axis_num (int, optional): The number of the axis that corresponds to the time axis in the
            displacement data array. Defaults to 1.

    Returns:
        np.array: A lateral position by time array containing the space-time plot, with reverb frames removed.
        np.array: A temporal positioning vector, adjusted to account for the removal of the reverb frames.
    """
    
    if not num_refs is None:
        
        reverb_idx = num_refs + 1
        
    else: # The number of reference frames is not provided, need to search for the true reference frame (all 0s)
        
        ref_frame = np.argmin(np.sum(np.abs(disp_plane), axis=lat_axis_num))
        reverb_idx = ref_frame + 2
        
    disp_plane_slicer_after_reverb = [slice(disp_plane.shape[ind]) for ind in range(len(disp_plane.shape))]
    disp_plane_slicer_after_reverb[t_axis_num] = slice(reverb_idx+n_additional_steps_to_remove, None)
    
    disp_plane = disp_plane[tuple(disp_plane_slicer_after_reverb)]
    
    t_ms = t_ms[reverb_idx+n_additional_steps_to_remove:]
    
    # Remove displacement offset at time 0 from all future frame
    
    disp_slicer_first_frame = [slice(disp_plane.shape[ind]) for ind in range(len(disp_plane.shape))]
    disp_slicer_first_frame[t_axis_num] = slice(0,1)
    
    disp_plane = disp_plane - disp_plane[tuple(disp_slicer_first_frame)]
    
    return disp_plane, t_ms

def low_pass_filter_plane_in_time_2D(disp_plane:np.array, t_ms:np.array, LPF_CO_kHz:float,
                                     filter_order:int=3, t_ax_num:int=1):
    """ Applies a 2D butter filter to the space-time plot.

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        t_ms (np.array): The vector containing the temporal information corresponding to the space-time plot.
        LPF_CO_kHz (float): The low-pass cutoff for the applied filter, in kHz.
        filter_order (int, optional): The filter order. Defaults to 3.
        t_ax_num (int, optional): The number of the axis that corresponds to the time axis in the
            displacement data array. Defaults to 1.

    Returns:
        np.array: The space-time plot, filtered by the butter filter.
    """
    
    dt_ms = np.mean(np.diff(t_ms))
    Fs_kz = 1/dt_ms
    
    Wn = (LPF_CO_kHz/(Fs_kz/2))
    
    b,a = sps.butter(filter_order, Wn)
    
    return sps.filtfilt(b, a, disp_plane, axis=t_ax_num)

def add_zeros_before_push_frame(disp_plane:np.array, t_ms:np.array, min_time_ms:float,
                                t_ax_num:int=1):
    """ Adds zeros to all time points prior to the push frame in the space-time plot.

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        t_ms (np.array): The vector containing the temporal information corresponding to the space-time plot.
        min_time_ms (float): The minimum time point for which zeros will be appended to the space-time plot
            to reach the minimum time point.
        t_ax_num (int, optional): The number of the axis that corresponds to the time axis in the
            displacement data array. Defaults to 1.
            
    Raises:
        ValueError: Raises a Value Error in the event that the minimum time specified to reach is larger than the
            current minimum time in the time array, indicating that zeros do not need to be appended to the
            space-time plot. 

    Returns:
        np.array: The space-time plot, with zeros appended to achieve a minimum time value prior to the
            push event.
        np.array: The new temporal vector, adjusted to account for the added zero's frames. 
    """
    
    dt_ms = np.mean(np.diff(t_ms))
    n_frames = np.floor((t_ms[0]-min_time_ms)/dt_ms)
    
    if n_frames < 0:
        raise ValueError("Something has gone wrong here.")
    
    t_ms = np.arange(min_time_ms, t_ms[-1]+dt_ms/2, dt_ms)
    
    zeros_shape = list(disp_plane.shape)
    zeros_shape[t_ax_num] = int(n_frames)
    
    zeros_lines = np.zeros(zeros_shape)
    
    return np.concatenate((zeros_lines, disp_plane), axis=t_ax_num), t_ms

def upsample_disp_plane_over_axis(disp_plane:np.array, current_ax:np.array, desired_PRF:float,
                                  upsamp_ax_num:int=1):
    """ Upsamples the given space-time plot over the desired axis, using a cubic spline interpolator.

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        current_ax (np.array): The vector containing positioning information (spacial or temporal)
            corresponding to the axis over which the space-time plot is going to be upsampled.
        desired_PRF (float): The desired final PRF of the upsampled data. 
        upsamp_ax_num (int, optional): The axis over which the plane should be upsampled. Defaults to 1,
            corresponding to the assumed temporal axis of the space-time plot.
            
    Raises:
        ValueError: Raises a Value Error in the event that the axis vector over which to upsample the
            plane is not the same length as the dimensions of the array in the data plane.

    Returns:
        np.array: The plane, upsampled to the desired PRF.
        np.array: A vector containing positioning information for the upsampled axis.
    """
    
    if disp_plane.shape[upsamp_ax_num] != current_ax.shape[0]:
        raise ValueError("The length of the unit axis and the axis over which to"+
                         " upsample the disp plane are not the same. Please"+
                         " double check and correct your inputs.")
    
    interpolator = spi.CubicSpline(current_ax, disp_plane, axis=upsamp_ax_num)
    
    experimental_PRF = np.mean(np.diff(current_ax))
    upsampling_factor = round(desired_PRF/experimental_PRF)
    
    new_ax_spacing = (1/experimental_PRF)/upsampling_factor
    upsampl_ax = np.arange(current_ax[0], current_ax[-1]+new_ax_spacing, new_ax_spacing)
    
    return interpolator(upsampl_ax), upsampl_ax

def apply_tukey_window_over_ax_edges(disp_plane:np.array, ax:np.array, tukey_alpha:float,
                                     ax_num:int=1):
    """ Applies a Tukey window over the specified plane axis.

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        ax (np.array): The vector containing positioning information (spacial or temporal)
            corresponding to the axis over which the Tukey window is going to be applied.
        tukey_alpha (float): The alpha value for the Tukey window.
        ax_num (int, optional): The axis over which the Tukey window should be applied. 
            Defaults to 1.
            
    Raises:
        ValueError: Raises a Value Error in the event that the axis vector over which to apply the
            Tukey window to the plane is not the same length as the dimensions of the array in
            the data plane.

    Returns:
        np.array: The space-time plot with the applied Tukey window.
    """
    
    if disp_plane.shape[ax_num] != ax.shape[0]:
        raise ValueError("The length of the unit axis and the axis over which to"+
                         " apply the Tukey window the disp plane are not the same. Please"+
                         " double check and correct your inputs.")
    
    mask = sps.windows.tukey(ax.shape[0], tukey_alpha)
    
    if ax_num == 0: 
        # Transposes the mask if the Tukey window is intended to be applied in the 0th axis.
        mask = mask[:, np.newaxis].T
        
    return disp_plane * mask

def differentiate_disp_plane_to_generate_vel_plane(disp_plane:np.array, t_ax:np.array, t_ax_num:int=1):
    """ Differentiates the space-time plot over the temporal axis in order to create a velocity-data 
        space-time plot.

    Args:
        disp_plane (np.array): A numpy array containing a space-time plot.
        t_ax (np.array): The vector containing the temporal information corresponding to the space-time plot.
        t_ax_num (int, optional): The number of the axis position containing temporal information
            in the spacetime plot. Defaults to 1.

    Returns:
        np.array: The velocity-data space-time plot.
        np.array: The vector containing temporal information for the velocity-data space-time plot.
    """
    
    vel_plane = np.diff(disp_plane, axis=t_ax_num)
    vel_t = (t_ax[:-1] + t_ax[1:])/2
    
    return vel_plane, vel_t