"""

DataAcqPar.py

This script contains the class definition and helper functions for QIBA DataAcqPar Objects.
These objects hold the acquisition parameters for the pipeline, and are necessary for its performance.

Author: Kaden Bock (03/29/2025)

"""

import numpy as np
import pymatreader
import os
import struct

class DataAcqPar:
    def __init__(self, directory_path:str, filename:str, push_V:float, push_freq:float, num_push_cycles:int, 
                 push_focus_mm:float, speed_of_sound:float, track_freq:float, NumEns:int, NumRefs:int, 
                 lastBmodeTransmit:int, T:np.array, T_indx:np.array, track_lat_axis:np.array,
                 track_axial_axis:np.array, pixels_per_track_frame:int):
        """ Initialization function to generate DataAcqPar objects, which store information about the
            acquisition parameters used to gather the SWEI data.

        Args:
            directory_path (str): The directory where the acquisition is stored.
            filename (str): The specific filename of the currently analyzed acquisition.
            push_V (float): The push voltage used to generate the acoustic radiation force.
            push_freq (float): The frequency of the acoustic radiation force push.
            num_push_cycles (int): The number of cycles for the push event.
            push_focus_mm (float): The axial focal depth of the push event.
            speed_of_sound (float): The assumed speed of sound for the system and RoI.
            track_freq (float): The frequency used for the tracking images.
            NumEns (int): The number of tracking frames gathered.
            NumRefs (int): The number of reference frames (prior to the push) gathered.
            lastBmodeTransmit (int): The number of the last Bmode Transmit in the Tx struct.
            T (np.array): The temporal axis corresponding to when tracking frames are gathered.
            T_indx (np.array): The indexes of the temporal axis.
            track_lat_axis (np.array): The tracking frames' lateral axis.
            track_axial_axis (np.array): The tracking frames' axial axis.
            pixels_per_track_frame (int): The total number pixels contained within a single
                tracking frame.
        """
        
        self.directory_path = directory_path
        self.full_param_file_name = filename
        self.filestamp = filename[:14]
        self.push_V = push_V
        self.push_freq_MHz = push_freq
        self.num_push_cycles = num_push_cycles
        self.push_focus_mm = push_focus_mm
        self.c_mps = speed_of_sound
        self.track_freq_MHz = track_freq
        self.num_ens = int(NumEns)
        self.num_refs = NumRefs
        self.lastBmodeTransmit = lastBmodeTransmit
        self.T = T
        self.T_idx = T_indx
        
        self.lat_axis_mm = track_lat_axis
        self.ax_axis_mm = track_axial_axis
        
        self.dimensions = (self.ax_axis_mm.size, self.lat_axis_mm.size, self.num_ens)
        
        self.pixels_per_track_frame = pixels_per_track_frame
        return
    
    def load_IQData(self):
        """ Loads the IQ Data associated with a given acquisition, based on the given
            acquisition parameters.

        Returns:
            np.array: The in-phase component of the IQ-Data.
            np.array: The quadrature component of the IQ-Data.
        """
        
        file_root = os.path.join(self.directory_path, self.filestamp)
        
        num_vals_to_read = int(self.pixels_per_track_frame*self.num_ens)
        
        IData = np.fromfile(file_root+'_IQreal.bin', dtype=np.int32, count=num_vals_to_read)
        QData = np.fromfile(file_root+'_IQimag.bin', dtype=np.int32, count=num_vals_to_read)
        
        IData = np.reshape(IData, self.dimensions, 'F')
        QData = np.reshape(QData, self.dimensions, 'F')
        
        return IData, QData
    
    def calc_num_pixels_for_num_WL(self, num_WL:float):
        
        WL_mm = (self.c_mps/(self.track_freq_MHz*1e6))*1E3 # mm
        
        mm_per_pixel = self.ax_axis_mm[1] - self.ax_axis_mm[0]
        
        return int((WL_mm*num_WL)/mm_per_pixel)
    
    
def create_DataAcqPar_from_QIBA_file(filepath:str):
    """Creates a DataAcqPar Object from the specified filepath.
    
    Currently, only .mat files generated from QIBA acquisition scripts
    are supported. 

    Args:
        filepath (str): A string specifying the location of the parameter file

    Returns:
        DataAcqPar: A DataAcqPar object which holds the information stored in the
        associated parameter file. 
    """
    
    directory_path, filestamp = check_valid_param_filepath(filepath)
    
    output_dict = load_acq_params_from_mat_file(filepath)
    
    push_V = output_dict['Vpush']
    PData = output_dict['PData']
    TW = output_dict['TW']
    Resource = output_dict['Resource']
    NumEns = int(output_dict['ne'])
    NumRefs = int(output_dict['nrefs'])
    lastBmodeTransmit = output_dict['lastBmodeTransmit']
    T = output_dict['T']
    T_idx = output_dict['T_idx']
    
    speed_of_sound = Resource['Parameters']['speedOfSound']
    
    track_freq = TW['Parameters'][0][0]
    
    push_freq = TW['Parameters'][1][0]
    num_push_cycles = TW['Parameters'][1][2]/2
    push_focus_wl = output_dict['TX']['focus'][int(output_dict['lastBmodeTransmit']+1)]
    push_focus_mm = push_focus_wl*(speed_of_sound/(track_freq*1E6))*1000
    
    track_lat_ax, track_axial_ax = parse_PData_for_tracking_axis_QIBA(PData, speed_of_sound, track_freq)
    
    pixels_per_track_frame = PData['Size'][1][0]*PData['Size'][1][1]

    return DataAcqPar(directory_path, filestamp, push_V, push_freq, num_push_cycles, push_focus_mm, speed_of_sound, track_freq, 
                     NumEns, NumRefs, lastBmodeTransmit, T, T_idx, track_lat_ax, track_axial_ax, pixels_per_track_frame)
        
        
def check_valid_param_filepath(filepath):
    """Checks input parameter filepath strings to ensure they exist and are valid
    parameter file types.
    
    Currently, only .mat files are supported. If a filepath extension is missing,
    it is assumed that the intended file is a .mat file, and the .mat extension is
    added.

    Args:
        filepath (str): A string specifying the location of the parameter file
        
    Raises:
        TypeError: In the event that an unsupported file type is specified.
        FileNotFoundError: In the event that the specified file cannot be found.
        
    Returns:
        Str: The filepath, with an added extension if needed.
    """
    
    file_ex = filepath.split(".")[-1]
    
    if file_ex != "mat":
        if "." not in filepath:
            Warning("No file extension detected, attempting to check if .mat file exists!")
            filepath = filepath + ".mat"
        else:
            raise TypeError("Incorrect file type specified, at this time only .mat \
                files can be correctly parsed!")
        
    if not os.path.isfile(filepath):
        raise FileNotFoundError("The specified acquisition parameters file was \
                                not able to be found, please check the input filepath.")
                    
    filepath = os.path.normpath(filepath)
    
    directory_path, filestamp = os.path.split(filepath)
    
    return directory_path, filestamp

def load_acq_params_from_mat_file(filepath):
    """Loads the necessary DataAcqPar params from the associated filepath.
    
    This function specifically loads .mat parameter files, and checks to
    ensure that the loaded parameters correctly contain all necessary
    parameters.

    Args:
        filepath (str): A string specifying the location of the parameter file

    Raises:
        ValueError: In the event that a necessary parameter cannot be found in
        the loaded parameter file.

    Returns:
        dict: A dictionary containing all loaded data from the parameter file.
    """
    expected_keys = ['T_idx', 'Vvalue', 'T', 'nrefs', 'Vpush', 'PData', 'TX', 'Trans', 
                     'ne', 'TW', 'Resource', 'lastBmodeEvent', 'Receive', 
                     'lastBmodeTransmit', 'pushAngleDegree', 'lastBmodeReceive']
    
    loaded_dict = pymatreader.read_mat(filepath)
    
    loaded_keys = loaded_dict.keys()
    
    for key in expected_keys:
        if key not in loaded_keys:
            raise ValueError("The required parameter {} was not found in the loaded \
                             parameters file.".format(key))
    
    return loaded_dict

def parse_PData_for_tracking_axis_QIBA(PData:dict, c:float, track_freq_MHz:float):
    """ Parses the PData struct in order to generate the lateral and axial axises for the tracking data.

    Args:
        PData (dict): The PData struct from Verasonics.
        c (float): The speed of sound used for reconstruction.
        track_freq_MHz (float): The tracking imaging frequency, in mHz.

    Returns:
        np.array: A numpy array containing the lateral positioning information for the tracking images.
        np.array: A numpy array containing the axial positioning information for the tracking images.
    """
    track_freq_Hz = track_freq_MHz*1e6
    lat_axis  = (np.arange(PData['Size'][1][1])*(PData['PDelta'][1][0])*(c/track_freq_Hz)+PData['Origin'][1][1])*1000
    lat_axis = lat_axis - np.average(lat_axis) # Shift the lat axis to center it at 0
    
    ax_axis = (np.arange(PData['Size'][1][0])*(PData['PDelta'][1][2])*(c/track_freq_Hz)+PData['Origin'][1][2])*1000
    
    return lat_axis, ax_axis
    
    