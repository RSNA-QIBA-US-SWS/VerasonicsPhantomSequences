
%==========================================================================

function AnalyzeAllAcquisitions

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
%   computeDIsplacementsSamTrack
%   phase_unwrap
%   AnalyzeARFIdata
%   SetupARIFdataPlane
%   Construct2DFTandCalcPhVel
%   MakePlot

    [filepath, name, ext] = fileparts(mfilename('fullpath'));
    addpath(filepath)
    
    paramFiles = dir('*_parameters.mat');       % get parameter files in directions

    for i = 1:length(paramFiles)    % for each file, get filestamp and compute displacements

        disp(['Processing ' num2str(i) ' of ' num2str(length(paramFiles))])     % comment to suppress output
        filestamp = paramFiles(i).name(1:14);
        genDispMTL(filestamp)
        AnalyzeARFIdata(filestamp);
    end

    MakePlot
end

%==========================================================================
