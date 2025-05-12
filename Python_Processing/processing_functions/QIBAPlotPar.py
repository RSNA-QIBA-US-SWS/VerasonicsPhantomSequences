"""

QIBAPlotPar.py

This script contains the class definition and helper functions for QIBA PlotPar Objects.
These objects hold the plotting parameters for the pipeline, and are necessary for plotting final outputs.

Author: Kaden Bock (04/20/2025)

"""

class QIBAPlotPar:
    def __init__(self, phantom_ID:str, max_plot_speed:float, max_plot_freq:float,
                 error_bar_frac:float, conf_int_factor:float, max_freq_for_phVel_plot:float,
                 low_freq_cutoff_pos_min:float, kr_thresh:float, max_phVel_speed:float,
                 th_factor:float, max_valid_gSWS:float, min_valid_gSWS:float):
        """ Initializes an object which holds information for plotting QIBA output plots.

        Args:
            phantom_ID (str): The name of the phantom being tested, for the plot title.
            max_plot_speed (float): The upper limit on the y-axis of the output plots.
            max_plot_freq (float): The upper limit on the x-axis of the output plots.
            error_bar_frac (float): The fraction of the mean that should be used for error bars.
            conf_int_factor (float): The factor used to calculate error bars.
            max_freq_for_phVel_plot (float): The maximum frequency for the phase velocity plots.
            low_freq_cutoff_pos_min (float): The low frequency cut-off for phase velocity displays.
            kr_thresh (float): The multiplicative threshold for the kr values.
            max_phVel_speed (float): The maximum speed for the phase velocity plot.
            th_factor (float): The factor used to determine valid phase velocities for plotting.
            max_valid_gSWS (float): The maximum value allowed for "valid" group shear wave speeds.
            min_valid_gSWS (float): The minimum value allowed for "valid" group shear wave speeds.
        """
        
        self.phantom_ID = phantom_ID
        self.max_plot_speed = max_plot_speed
        self.max_plot_freq = max_plot_freq
        self.error_bar_frac = error_bar_frac
        self.conf_int_factor = conf_int_factor
        self.max_freq_for_phVel_plot = max_freq_for_phVel_plot
        self.low_freq_cutoff_pos_min = low_freq_cutoff_pos_min
        self.kr_thresh = kr_thresh
        self.max_phVel_speed = max_phVel_speed
        self.th_factor = th_factor
        self.max_valid_gSWS = max_valid_gSWS
        self.min_valid_gSWS = min_valid_gSWS
        
        pass
    