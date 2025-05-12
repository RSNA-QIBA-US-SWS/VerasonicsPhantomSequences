"""

AnalysisPar.py

This script contains the class definition and helper functions for QIBA AnalysisPar Objects.
These objects hold the analysis parameters for the pipeline, and are necessary for its performance.

Author: Kaden Bock (03/29/2025)

"""

import numpy as np

class AnalysisPar:
    def __init__(self, save_dir:str, kasai_kernel_WL:float, ax_proc_depth:float, 
                 DOF_mm:float, min_lat_mm:float, max_lat_mm:float, min_time_ms:float, 
                 max_time_ms:float, rolloff_time_ms:float, n_steps_to_rm:int, 
                 LPF_co_kHz:float, des_PRF_kHz:float, min_PV_freq_hz:float,
                 max_PV_freq_hz:float, PV_analysis_step:float):
        """ Initialization function for objects that store QIBA analysis parameters.

        Args:
            save_dir (str): The directory in which generated outputs should be stored.
            kasai_kernel_WL (float): The kasai axial kernel length, specified as a multiple
                of the imaging wavelength.
            ax_proc_depth (float): The axial depth at which the space-time plot should
                be created. Should be in the same units as the axial positioning axis.
            DOF_mm (float): The axial depth of field over which the axial information
                should be averaged.
            min_lat_mm (float): The min lateral position (mm) to start analysis.
            max_lat_mm (_type_): The max lateral position (mm) to end analysis.
            min_time_ms (_type_): The minimum time (ms) to zero-pad to before the push.
            max_time_ms (_type_): The maximum time (ms) to allow over the entire space-time plot.
            rolloff_time_ms (_type_): The time (ms) for rolloff at early and late times in the space-time plot.
            n_steps_to_rm (_type_): The number of reverb steps to remove after the push.
            LPF_co_kHz (_type_): The temporal low-pass filter cutoff frequency in kHz
            des_PRF_kHz (_type_): The desired temporal PRF after uspampling.
            min_PV_freq_hz (_type_): The minimum phase velocity to extract.
            max_PV_freq_hz (_type_): The maximum phase velocity to extract. 
            PV_analysis_step (_type_): The step size between extracted phase velocities.
        """
        
        self.save_dir = save_dir
        
        self.kasai_kernel_size_WL = kasai_kernel_WL
    
        self.axial_processing_depth = ax_proc_depth
        self.Depth_of_Field_to_average_mm = DOF_mm # depth of field to average (mm)
            
        self.min_lat_mm         = min_lat_mm          # min lateral position (mm) to start analysis
        self.max_lat_mm         = max_lat_mm       # max lateral position (mm) to end analysis
        self.min_time_ms        = min_time_ms        # maximum time (ms) to use in analysis
        self.max_time_ms        = max_time_ms       # minimum time (ms) to zero-pad before push
        self.rolloff_time_ms    = rolloff_time_ms          # time (ms) for rolloff at early and late times
        self.n_steps_to_remove   = n_steps_to_rm           # number of reverb steps to remove
        self.LPF_cutoff_kHz     = LPF_co_kHz          # low-pass filter cutoff frequency (kHz)
        self.desired_PRF_kHz    = des_PRF_kHz          # desired PRF (kHz) after upsampling in time 
        self.min_phase_vel_freq_Hz = min_PV_freq_hz   # frequencies (Hz) for phase velocity measurements
        self.max_phase_vel_freq_Hz = max_PV_freq_hz
        self.phase_vel_analysis_step_Hz = PV_analysis_step
        
        total_time_width = max_time_ms - min_time_ms
        tukey_roll_off_width = 2*rolloff_time_ms
        tukey_alpha = tukey_roll_off_width/total_time_width
        
        if tukey_alpha > 1:
            print("WARNING: The requested roll_off_time_ms is larger than the" +
                          " total time window width. The Tukey window has been adjusted"
                          + " to one full cosine period.")
            
        self.tukey_alpha = tukey_alpha
        
        self.freqs_to_analyze = self.generate_freqs_to_analyze_vec()
        
        self.acq_params = None
        pass
    
    def generate_freqs_to_analyze_vec(self):
        """ Generates the frequency vector over which phase velocities will be
            be extracted.

        Returns:
            np.array: The frequency vector over which phase velocities will be
                be extracted.
        """
        
        start = self.min_phase_vel_freq_Hz
        stop = self.max_phase_vel_freq_Hz
        dfreq = self.phase_vel_analysis_step_Hz
        
        return np.arange(start, stop+dfreq/2, dfreq)
    
def generate_AnalysisPar_obj_from_input_dict(analysis_par_dict:dict):
    """ Generates an AnalysisPar object that holds all of the analysis parameters from
        an input dictionary containing all the necessary information.
        
    Args:
        analysis_par_dict (dict): A dictionary that contains all of the necessary analysis
            parameters.

    Raises:
        KeyError: In the event that a necessary analysis parameter is missing from the
            input dictionary. 

    Returns:
        AnalysisPar: A QIBA.py AnalysisPar object that holds all the necessary parameters
            for QIBA processing.
    """
    
    expected_keys = ['kasai_kernel_size_WL', 'save_dir', 'axial_processing_depth',
                     'Depth_of_Field_to_average_mm', 'min_lat_mm', 'max_lat_mm', 'max_time_ms',
                     'min_time_ms','rolloff_time_ms','n_steps_to_remove','LPF_cutoff_kHz', 'desired_PRF_kHz',
                     'min_phase_vel_freq_Hz','max_phase_vel_freq_Hz','phase_vel_analysis_step_Hz']
    
    actual_keys = analysis_par_dict.keys()
    
    for key in expected_keys:
        try:
            assert(key in actual_keys)
        except AssertionError:
            raise KeyError("Please ensure that the processing parameter {} is defined.".format(key))
        
    ks_kn_wl = analysis_par_dict['kasai_kernel_size_WL']
    save_dir = analysis_par_dict['save_dir']
    proc_depth = analysis_par_dict['axial_processing_depth']
    DOF_mm = analysis_par_dict['Depth_of_Field_to_average_mm']
    min_lat_mm = analysis_par_dict['min_lat_mm']
    max_lat_mm = analysis_par_dict['max_lat_mm']
    max_tms = analysis_par_dict['max_time_ms']
    min_tms = analysis_par_dict['min_time_ms']
    ro_tms = analysis_par_dict['rolloff_time_ms']
    n_steps = analysis_par_dict['n_steps_to_remove']
    LPF_cutoff = analysis_par_dict['LPF_cutoff_kHz']
    desired_PRF = analysis_par_dict['desired_PRF_kHz']
    min_ph = analysis_par_dict['min_phase_vel_freq_Hz']
    max_ph = analysis_par_dict['max_phase_vel_freq_Hz']
    ph_vel_step = analysis_par_dict['phase_vel_analysis_step_Hz']
    
    par_obj = AnalysisPar(save_dir, ks_kn_wl, proc_depth, DOF_mm, min_lat_mm,
                          max_lat_mm, min_tms, max_tms, ro_tms, n_steps,
                          LPF_cutoff, desired_PRF, min_ph, max_ph, ph_vel_step)
    return par_obj