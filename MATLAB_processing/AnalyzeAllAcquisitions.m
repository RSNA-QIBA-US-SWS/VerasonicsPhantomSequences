
%==========================================================================

function AnalyzeAllAcquisitions

    clc;close all;
    addpath(pwd)

    % Run this script from a directory containing the raw data to compute
    % displacements for all of the timestamped files in that directory.
    % 
    % Timestep should have 14 characters (YYYYMMDDHHMMSS), e.g. 20171015104934.
    % For each timestep, you need three files:
    %   - [timestepNumber]_IQreal.bin
    %   - [timestepNumber]_IQimag.bin
    %   - [timestepNumber]_parameters.mat
    %
    % The following scripts must be in the same directory, or use addpath
    %   genDispMTL, 
    %   kasai_algorithm,
    %   phase_unwrap
    %   AnalyzeARFIdata
    %   SetupARIFdataPlane
    %   Construct2DFTandCalcPhVel
    %   MakePlotAndSaveCSVfile

    !dataDir = 'data_directory';    % CHANGE ME to a directory string or "pwd" for the current directory
    cd(dataDir)

       % analysis directory definition

    par.analysisDir      = dataDir;         % directory of verasonics data, and output plot save directory

        % parameters for processing arfidata

    par.DOFmm            = 2;           % depth of field to average (mm)
    par.minLatMM         = 6;           % min lateral position (mm) to start analysis
    par.maxLatMM         = 20;          % max lateral position (mm) to end analysis
    par.maxTimeMS        = 40;          % maximum time (ms) to use in analysis
    par.minTimeMS        = -3;          % minimum time (ms) to zero-pad before push
    par.rolloffTimeMS    = 25;          % time (ms) for rolloff at early and late times
    par.nStepsToRemove   = 2;           % number of reverb steps to remove
    par.LPFcutoffKHz     = 1;           % low-pass filter cutoff frequency (kHz)
    par.desiredPRFkHz    = 20;          % desired PRF (kHz) after upsampling in time 
    par.freqsToAnalyzeHz = 20:10:800;   % frequencies (Hz) for phase velocity measurements
    par.plot_int_fig     = 0;           % flag for plotting intermediate output figures (shearwave propagation animation and spacetime plot)
    par.swsest           = 'ttp';

        % parameters for plotting and saving output
    
    par.phantomID        = 'phantom_id';  % phantom ID string
    par.maxPlotSpeed     = 6;         % maximum speed for gSWS and phVel plots
    par.maxplotfreq      = 500;         % max frequency (Hz) for plot
    par.fracErrorBar     = 0.3;         % fixed fraction for error bars
    par.CIfactor         = 1.96;        % confidence interval for error bars
    par.fmax             = 600;         % maximum frequency (Hz) to plot
    par.rmin             = 0.004;       % minimum lateral position (m) for low freq cutoff
    par.krThreshold      = 1.5;         % kr threshold for low freq cutoff
    par.maxSpeed         = 8;           % max speed for "good" result
    par.thFactor         = 2;           % phVel accept range = +/- thfactor * std
    par.maxValidGSWS     = 6;           % maximum group speed where results are valid
    par.minValidGSWS     = 0.5;         % minimum group speed where results are valid

   % [filepath, name, ext] = fileparts(mfilename('fullpath'));
   % addpath(filepath)
   
    paramFiles = dir([dataDir '/*_parameters.mat']);    % get parameter files

    for i = 1:length(paramFiles)    % for each file, get filestamp and compute displacements
        disp(['Processing ' num2str(i) ' of ' num2str(length(paramFiles))])     % comment to suppress output
        filestamp = paramFiles(i).name(1:14);
        par.filestamp = filestamp;
        genDispMTL(filestamp, par)
        AnalyzeARFIdata(filestamp,par);
    end

    errorFlag = MakePlotAndSaveOutputCSVfile(par);

end

%==========================================================================
