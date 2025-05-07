"""

genDisp.py

This file contains all functions needed to generate displacement estimates from IQ Data files.

Author: Kaden Bock (03/27/2025)

"""

import numpy as np
import scipy.signal as sps

def calc_particle_displacement_using_kasai(IQData: np.array, num_refs: int, kernel_size_pixels: int,
                                           c_mps: float, track_freq_Mhz: float, reference_mode: str = 'fixed',
                                           kasai_axis: int = 0, big_time_axis: int = 2):
    """ Calculates particle motion using Kasai's algorithm on demodulated IQ data.

        For more information on Kasai's algorithm, refer to: https://doi.org/10.1109/tuffc.2006.1642509

    Args:
        IQData (np.array): IQData array containing the complex beamformed tracking frame data. 
            * It is assumed that the slow time axis is the 3rd axis in this array.
        num_refs (int): The number of reference frames taken before the push.
            * If the number of reference frames is greater than 1, then the additional reference 
                frames are tracked in reverse.
        kernel_length (int): The length of the kasai kernel, in number of pixels
        c (float): The assumed speed of sound of material, in meters/second
        track_freq (float): The frequency of the tracking wave, in MHz.
        reference_mode (str, optional): Allows the user to select between using a progressive
            or fixed reference scheme. Defaults to 'fixed', with the other implemented option
            being "progressive."
        kasai_axis (int, optional): Selects which axis of the IQData array to calculate phase shifts
            over. This should be set to the same axis that axial indexing is stored in. 
            Defaults to 0.
        big_time_axis (int, optional): Selects which axis of the IQData array is the big time-axis.
            This should be set to the same axis that temporal frame indexing is stored in. 
            Defaults to 2.

    Raises:
        IndexError: Raises an index error if no reference frame prior to the push is provided.
        NotImplementedError: Raises NotImplementedErrors in the event that the the IQData
            matrix is not a third-rank matrix, or if the reference mode specified is not
            either "progressive" or "fixed"

    Returns:
        np.array: A numpy array containing the calculated particle motion for each pixel through the
            indicated kasai axis.
    """

    if num_refs < 1:
        raise IndexError("There needs to be at least one track frame prior to the push in order" +
                         " to be able to track particle motion after the push")

    if len(IQData.shape) != 3:
        raise NotImplementedError(
            "Analyzing IQData on not-third order input matrices is not supported.")

    # Create a list of slicers to slice through the IQ Data array
    IQData_slicer = [slice(IQData.shape[ind])
                     for ind in range(len(IQData.shape))]

    ref_idx = num_refs - 1  # Convert to python 0-indexing

    if reference_mode == 'progressive':
        # ref_data = IQData[:,:,ref_idx:-1]
        IQData_slicer[big_time_axis] = slice(ref_idx, -1)
        # Store a copy of the reference frame
        ref_data = IQData[tuple(IQData_slicer)]
    elif reference_mode == 'fixed':
        # ref_data = IQData[:,:,ref_idx:num_refs]
        IQData_slicer[big_time_axis] = slice(ref_idx, num_refs)
        # Slice the array into progressive reference frames
        ref_data = IQData[tuple(IQData_slicer)]
    else:
        raise NotImplementedError(
            "The requested reference mode is not supported by this motion estimator.")

    # kData = IQData[:,:,num_refs:]
    IQData_slicer[big_time_axis] = slice(num_refs, None)
    kData = IQData[tuple(IQData_slicer)]

    kasai_numerator = (kData.imag*ref_data.real) - (kData.real*ref_data.imag)
    kasai_denominator = (kData.real*ref_data.real) + (kData.imag*ref_data.imag)

    # Create the convolution kernel
    kernel_shape = np.ones(len(kData.shape), dtype=int)
    kernel_shape[kasai_axis] = kernel_size_pixels
    kernel = np.ones(kernel_shape)/kernel_size_pixels

    kasai_numerator = sps.fftconvolve(
        kasai_numerator, kernel, mode='same', axes=kasai_axis)
    kasai_denominator = sps.fftconvolve(
        kasai_denominator, kernel, mode='same', axes=kasai_axis)

    disp_est = np.arctan2(kasai_numerator, kasai_denominator) * \
        (c_mps/(4*track_freq_Mhz*np.pi))

    if num_refs > 1:
        # If there is more than one reference frame, track the additional reference frames in reverse
        # relative to the reference frame right before the push

        IQData_slicer[big_time_axis] = slice(0, num_refs)

        # additional_reference_frames = IQData[:,:, 0:num_refs]
        additional_reference_frames = IQData[tuple(IQData_slicer)]
        additional_reference_frames = np.flip(
            additional_reference_frames, axis=big_time_axis)

        pretrack_disp_est = -1 * calc_particle_displacement_using_kasai(additional_reference_frames, 1,
                                                                        kernel_size_pixels, c_mps,
                                                                        track_freq_Mhz, 'fixed', kasai_axis)

        # Re-reversing the frames to put them back in the correct order
        pretrack_disp_est = np.flip(pretrack_disp_est, axis=big_time_axis)

        disp_est = np.concatenate(
            (pretrack_disp_est, disp_est), axis=big_time_axis)

    else:
        # Add a frame of all 0s as a stand-in for the reference frame
        # (Since all motion is relative to the reference, the reference compared to itself has 0 motion)

        zeros_shape = list(IQData.shape)
        zeros_shape[big_time_axis] = 1
        ref_zero_frame = np.zeros(zeros_shape)

        disp_est = np.concatenate(
            (ref_zero_frame, disp_est), axis=big_time_axis)

    if reference_mode == "progressive":
        # If the motion estimates were generated using a progressive reference scheme, integrate the
        # estimates to get particle displacement

        disp_est = np.cumsum(disp_est, axis=big_time_axis)

    # Invert the calculated displacement such that it aligns with our imaging orientation
    return -1 * disp_est


def phase_unwrap_displacement_estimate(disp_est: np.array, c_mps: float, track_freq_MHz: float,
                                       first_frame_to_check: int, big_time_axis: int = 2,
                                       scale_factor: float = 1.0, total_iters: int = 1):
    """ Unwraps phase estimates in the event that the displacement estimate exceeds
        the imaging wavelength.

    Args:
        disp_est (np.array): A numpy containing the calculated particle motion for each pixel through the
            indicated kasai axis.
        c_mps (float): Assumed speed of sound of the system, in meters/sec
        track_freq_MHz (float): The wave frequency used for tracking images, in MHz
        first_frame_to_check (int): First big-time frame for which displacement estimates should be unwrapped
        big_time_axis (int, optional): Axis of the displacement estimates over which there are different tracking 
            frames. Defaults to 2.
        scale_factor (float, optional): Scaling factor by which the critical distance should be scaled. 
            Defaults to 1.0.
        total_iters (int, optional): Total number of times to unwrap phase estimates. Defaults to 1.

    Returns:
        np.array: The array of phase-unwrapped displacement estimations
    """

    crit_dist = (c_mps / (track_freq_MHz*4)) * scale_factor  # pi's cancel out

    disp_est_slicer = [slice(disp_est.shape[ind])
                       for ind in range(len(disp_est.shape))]

    disp_est_slicer[big_time_axis] = slice(0, first_frame_to_check+1)

    for iter in range(total_iters):
        diff = np.diff(disp_est, axis=big_time_axis, prepend=0)
        diff = -1*diff*(diff > crit_dist) + diff*(diff < -1*crit_dist)

        diff[tuple(disp_est_slicer)] = 0

        integ = np.cumsum(diff, axis=big_time_axis)*crit_dist

        disp_est = disp_est + integ*2

    return disp_est
