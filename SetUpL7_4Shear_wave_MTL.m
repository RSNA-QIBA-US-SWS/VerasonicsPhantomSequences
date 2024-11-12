clear all
global filedir
scriptName='SETUPL7_4Shear_wave_MTL';

%% Commonly changing variables
%filedir = '/home/verasonics/cloud/Vantage-4.2.0-2001220500';
filedir = 'C:\Users\Administrator\Documents\VantageNXT-1.1.13\';
cd(filedir)
saveChannelData = 0;
saveIQData = 1;
push_cycle  = 900;
push_focus = 25; % in mm
push_Fnum = 1.5;
npush = 1;
ne = 50;
nrefs = 5;
pushAngleDegree = 0;
%% Define basic parameters
m = 128; % Bmode lines
getPower = 0; % DO NOT delete, used by save_swei_data

Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;

P.numRays = m;
P.startDepth = 0;   % Acquisition depth in wavelengths
P.endDepth = 240;   % This should preferrably be a multiple of 128 samples.

%% Specify Trans structure array.
Trans.name = 'L7-4';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans(Trans);  % C5-2 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 76;  % set maximum high voltage limit for pulser supply.
TPC(5).maxHighVoltage = 76;
w = Resource.Parameters.speedOfSound/Trans.frequency/1000; % wavelength in mm
pushElementNum = round((push_focus/push_Fnum)/(Trans.spacing*w)/2)*2;

%% Specify PData structure array (only used for Bmode recon)
PData(1).PDelta = [Trans.spacing, 0, 1.0];
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil((P.numRays*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Rectangle',...
                    'Position',[0,0,P.startDepth],...
                    'width',Trans.spacing,...
                    'height',P.endDepth-P.startDepth)),1,128);
% - set position of regions to correspond to beam spacing.
for i = 1:128
    PData(1).Region(i).Shape.Position(1) = (-63.5 + (i-1))*Trans.spacing;
end
PData(1).Region = computeRegions(PData(1));

PData(2).PDelta = [Trans.spacing/2, 0, 0.25];
PData(2).Size(1) = ceil((P.endDepth-P.startDepth)/PData(2).PDelta(3));
PData(2).Size(2) = ceil(((Trans.numelements-pushElementNum)*Trans.spacing)/PData(2).PDelta(1))+1;
PData(2).Size(3) = 1;      % single image page
PData(2).Origin = [-Trans.spacing*(Trans.numelements-pushElementNum)/2,0,P.startDepth] ; % x,y,z of upper lft crnr
% [PData(2).Region,PData(2).numRegions]=createRegions(PData(2));

%% Specify resource buffers
% RcvBuffer stores channel data. Buffer 1 stores Bmode data
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2400*m; % max 4096 per axial line
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 2; 

Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = ne*2400;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = npush; %%% Change this number to acquire multiple focal zones

Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = npush;  % one intermediate buffer needed.
Resource.InterBuffer(1).rowsPerFrame = PData(2).Size(1);
Resource.InterBuffer(1).colsPerFrame = PData(2).Size(2);
Resource.InterBuffer(1).pagesPerFrame = ne;

% Image buffer saves reconstructed intensity data for Bmode image
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).rowsPerFrame = PData(1).Size(1); % this is for maximum depth
Resource.ImageBuffer(1).colsPerFrame = PData(1).Size(2);
Resource.ImageBuffer(1).numFrames = 1;

% Set up Bmode display window
Resource.DisplayWindow(1).Title = 'Image Display';
Resource.DisplayWindow(1).pdelta = 0.4;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);

%% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];   % A, B, C, D

% - Push waveform.
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,1,push_cycle*2,1];  %

%% Specify TX structure array.  
% Specify nr TX structure arrays. Transmit centered on element n in the array for event n.
txFocus = round(25/w);  % Initial transmit focus.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', txFocus, ... 
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, m+1+npush);

% Determine TX aperture based on focal point and desired f number.
txFNum = 2;  % set to desired f-number value for transmit (range: 1.0 - 20)
txNumEl=round((txFocus/txFNum)/Trans.spacing/2); % no. of elements in 1/2 aperture.
if txNumEl > (Trans.numelements/2 - 1), txNumEl = floor(Trans.numelements/2 - 1); end   
% txNumEl is the number of elements to include on each side of the
% center element, for the specified focus and sensitivity cutoff.
% Thus the full transmit aperture will be 2*txNumEl + 1 elements.
%display('Number of elements in transmit aperture:');
%disp(2*txNumEl+1);
           
% - Set event specific TX attributes.
for n = 1:m   % 128 transmit events
    % Set transmit Origins to positions of elements.
    TX(n).Origin = [(-63.5 + (n-1)*Trans.spacing), 0.0, 0.0];
    % Set transmit Apodization so (1 + 2*TXnumel) transmitters are active.
    lft = n - txNumEl;
    if lft < 1, lft = 1; end;
    rt = n + txNumEl;
    if rt > Trans.numelements, rt = Trans.numelements; end;
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end
lastBmodeTransmit = n;

% Track
n = m+1;
TX(n).Origin = [0.0 0.0 0.0];
TX(n).focus = 0.0;
TX(n).Apod = ones(1,Resource.Parameters.numRcvChannels); %All channels
TX(n).Delay = computeTXDelays(TX(n));

% Push
onele_push = pushElementNum;
offele_push = Trans.numelements-onele_push;
for ipush = 1:npush
    n = m+1+ipush;
    TX(n).waveform = 2;
    TX(n).Origin = [0.0 0.0 0.0];
    TX(n).focus = round(push_focus/w);
    centerElement(ipush) = 64;
    if centerElement(ipush) < onele_push/2
        TX(n).Apod(1:onele_push) = 1;
    elseif Resource.Parameters.numRcvChannels-centerElement(ipush) < onele_push/2
        TX(n).Apod(onele_push+1:end) = 1;
    else
        TX(n).Apod(centerElement(ipush)-onele_push/2+1:centerElement(ipush)+onele_push/2) = 1;
    end
    TX(n).Delay = computeTXDelays(TX(n));
end

%% Specify Receive structure arrays. 
% Compute the maximum receive path length, using the law of cosines.
maxAcqLength = sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2) - P.startDepth;
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW',...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, Resource.RcvBuffer(1).numFrames*m+ne*npush); %%% Change this number to acquire multiple focal zones
%                         'InputFilter',repmat([0.0036 0.0127 0.0066 -0.0881 -0.2595 0.6494],Resource.Parameters.numRcvChannels,1),...
% - Set event specific Receive attributes for inital Bmode.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = m*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:m
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
    end
end

lastBmodeReceive = m*(i-1)+j;  

Receive(lastBmodeReceive+1).callMediaFunc = 0;
for j = 1:ne*npush %%% Change this number to acquire multiple focal zones
    Receive(lastBmodeReceive+j).Apod(:) = 1.0;
    Receive(lastBmodeReceive+j).bufnum = 2;
    Receive(lastBmodeReceive+j).framenum = ceil(j/ne);
    if mod(j,ne) ~=0
        Receive(lastBmodeReceive+j).acqNum = mod(j,ne);
    else
        Receive(lastBmodeReceive+j).acqNum = ne;
    end
    Receive(lastBmodeReceive+j).startDepth = 0;
end

%% Specify TGC Waveform structure.
TGC.CntrlPts = [400,550,650,710,770,830,890,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Recon structure arrays.
Recon(1) = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'newFrameTimeout',1000,...
               'rcvBufFrame',-1, ...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [1,-1], ...
               'RINums', [1:m]);

for i = 1:npush
    Recon(i+1) = struct('senscutoff', 0.6, ...
        'pdatanum', 2, ...
        'newFrameTimeout',2000,...
        'rcvBufFrame',i,...
        'IntBufDest', [1,i], ...
        'ImgBufDest', [0,0], ...
        'RINums',(m+1+ne*(i-1):(m+ne*i))');
end           
           
%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 0, ...  % replace data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, m+ne*npush);
% - Set specific ReconInfo attributes.
for i = 1:m
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end

% - ReconInfo for ARFI/SWEI.
k = m; % k keeps track of index of last ReconInfo defined
% We need ne*nacqs ReconInfo structures for IQ reconstructions.
ReconInfo((k+1):(k+ne*npush)) = repmat(struct('mode', 3, ... % IQ output 
                   'txnum', m+1, ...
                   'rcvnum', lastBmodeReceive+1, ...
                   'regionnum', 1), 1, ne*npush);

for j = 1:ne*npush  % For each row in the column
    ReconInfo(k+j).rcvnum = lastBmodeReceive+j;
    if mod(j,ne) ~=0
        ReconInfo(k+j).pagenum = mod(j,ne);
    else
        ReconInfo(k+j).pagenum = ne;
    end
    ReconInfo(k+j).txnum = m+1;
end
%% Specify Process structure array.
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1, ...   % number of buffer to process.
                         'framenum',-1, ...   % frame number in src buffer (-1 => lastFrame)
                         'pdatanum',1, ...    % number of PData structure (defines output figure).
                         'norm',1, ...        % normalization method(1 means fixed)
                         'pgain',20.0, ...            % pgain is image processing gain
                         'persistMethod','none', ...
                         'persistLevel',0, ...
                         'interp',1, ...      % method of interpolation (1=4pt interp)
                         'compression',0.5, ...      % X^0.5 normalized to output word size
                         'mappingMode','full', ...
                         'display',1, ...     % display image after processing
                         'displayWindow',1, ...
                         'compressMethod', 'power',...
                         'compressFactor', 40}; %compressMethod and compressFactor are both required NXT commands that were not previously here
         
Process(2).classname = 'External';
Process(2).method = 'save_channel_data';
Process(2).Parameters = {'srcbuffer','receive',...
                         'srcbufnum',2,...
                         'srcframenum',0,...
                         'dstbuffer','none'};
                  
Process(3).classname = 'External';
Process(3).method = 'save_IQ_data';
Process(3).Parameters = {'srcbuffer','inter',...
                         'srcbufnum',1,...
                         'srcframenum',0,... %1
                         'dstbuffer','none'};
                     
Process(4).classname = 'External';
Process(4).method = 'rt_matlab';
Process(4).Parameters = {'srcbuffer','inter',...
                         'srcbufnum',1,...
                         'srcframenum',0,... %1
                         'dstbuffer','none'};  

Process(5).classname = 'External';
Process(5).method = 'switch_power_initial';
Process(5).Parameters = {'srcbuffer','none',...
                         'srcbufnum',1,...
                         'srcframenum',0,...
                         'dstbuffer','none'};         

%% Specify SeqControl structure arrays.
% - Change to Profile 1
SeqControl(1).command = 'setTPCProfile';
SeqControl(1).condition = 'immediate';
SeqControl(1).argument = 1;
% - Noop to allow time for charging external cap.
SeqControl(2).command = 'noop';
SeqControl(2).argument = 500000; % wait 100 msec.
% - Set time between rays
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 250;
% - Set time between frames
SeqControl(4).command = 'timeToNextAcq';
SeqControl(4).argument = 100000;
% - Return to Matlab
SeqControl(5).command = 'returnToMatlab';
% - Jump back to 3.
SeqControl(6).command = 'jump';
SeqControl(6).argument = 3;

% - ARFI timing
% - Change to Profile 5 (high power)
SeqControl(7).command = 'setTPCProfile';
SeqControl(7).condition = 'immediate';
SeqControl(7).argument = 5;
% - time between tracks
SeqControl(8).command = 'timeToNextAcq';
SeqControl(8).argument = 200;
% - Trigger out
SeqControl(9).command = 'triggerOut';
% - time between pushes
SeqControl(10).command = 'timeToNextAcq';
SeqControl(10).argument = 200;
% - time between last Bmode tx and jump back 1st Bmode tx to include time
% for setting power level
SeqControl(11).command = 'timeToNextAcq';
SeqControl(11).argument = 500500;
% - Jump back to start.
SeqControl(12).command = 'jump';
SeqControl(12).argument = 1;

nsc = 13;

% Specify Event structure arrays.
n = 1;

TTNAS=zeros(ne*3,1); %keep track of timing (TimeToNextAcq's)
T_i=1;             %current index

%% Start of Events
Event(n).info = 'Switch to profile 1.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

Event(n).info = 'noop for charging ext. cap.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 2;
n = n+1;

% B-Mode Loop!
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:m                  % Acquire rays
        Event(n).info = 'Acquire ray line';
        Event(n).tx = j; 
        Event(n).rcv = m*(i-1)+j;   
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 3;
        n = n+1;
    end
    % Replace last event's SeqControl for inter-frame timeToNextAcq.
    Event(n-1).seqControl = 4;
    
    Event(n).info = 'Transfer frame to host.';
    Event(n).tx = 0;        % no TX
    Event(n).rcv = 0;       % no Rcv
    Event(n).recon = 0;     % no Recon
    Event(n).process = 0; 
    Event(n).seqControl = nsc; 
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc = nsc+1;
    n = n+1;

    Event(n).info = 'reconstruct'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1;      % reconstruction
    Event(n).process = 0;    % process
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'process (Display B-Mode) and return to Matlab'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % reconstruction
    Event(n).process = 1;    % process
    Event(n).seqControl = 0;
    if floor(i/2) == i/2     % Exit to Matlab every xth frame reconstructed 
        Event(n).seqControl = 5; %'returnToMatlab';
    end
    n = n+1;
end

Event(n).info = 'Jump back to third event (stay at current power)';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0; 
Event(n).seqControl = 6;
n = n+1;

lastBmodeEvent = n;

% Start of ARFI Events!
% Switch to TPC profile 5 (high power) and allow time for charging ext. cap.
Event(n).info = 'Switch to profile 5.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;
n = n+1;

Event(n).info = 'noop for charging ext. cap.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 2;
n = n+1;

% ARFI loop
for i = 1:npush 
    %push 1 ensemble
    for j = 1:nrefs
        Event(n).info = 'Acquire reference data';
        Event(n).tx = m+1;
        Event(n).rcv = lastBmodeReceive+(i-1)*ne+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 8;
        if T_i == 1 %if this is the first one, add it to zero
            TTNAS(T_i)=TTNAS(T_i)+SeqControl(Event(n).seqControl).argument; T_i=T_i+1;
        else        %otherise add new time to previous sum
            TTNAS(T_i)=TTNAS(T_i-1)+SeqControl(Event(n).seqControl).argument; T_i=T_i+1;
        end
        n = n+1;
    end
    
    Event(n).info = 'Push transmit';
    Event(n).tx = m+1+i;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = 10;
    TTNAS(T_i)=TTNAS(T_i-1)+SeqControl(Event(n).seqControl).argument;
    n = n+1;
    
    for j = nrefs+1:ne
        Event(n).info = 'Acquire tracking data';
        Event(n).tx = m+1;
        Event(n).rcv = lastBmodeReceive+(i-1)*ne+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 8;
        if j == nrefs+1
            TTNAS(T_i)=TTNAS(T_i)+SeqControl(Event(n).seqControl).argument; T_i=T_i+1;
        else
            TTNAS(T_i)=TTNAS(T_i-1)+SeqControl(Event(n).seqControl).argument; T_i=T_i+1;
        end
        n = n+1;
    end
    
    Event(n-1).seqControl = 11; % modify last detect acquisition's seqControl for frame interval
    
    Event(n).info = 'transfer data to Host';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no process
    Event(n).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 1+i;      % reconstruction
    Event(n).process = 0;    % process
    Event(n).seqControl = 2; % NOOP for .1sec
    n = n+1;
    
end

if saveChannelData
    Event(n).info = 'save channel data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 2;
    Event(n).seqControl = 2;
    n = n+1;
end

if saveIQData
    Event(n).info = 'save IQ data';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 2;
    n = n+1;
end

Event(n).info = 'Back to Matlab';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % no process
Event(n).seqControl = 5; % Back to Matlab
n = n+1;

Event(n).info = 'Jump back to live B-mode';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     
Event(n).process = 0;
Event(n).seqControl = 12;
n = n+1;

% Motion filter timing. time zero at the first push
T=(TTNAS(1:ne)-TTNAS(nrefs)-SeqControl(8).argument)/1000;
i=find(T>4.0,1,'first');
% i=length(T);
T_idx = [1:nrefs i':length(T)];

%% User specified UI Control Elements
% - Sensitivity Cutoff
sensx = 170;
sensy = 420;
UI(1).Control = {'Style','text',...        % popupmenu gives list of choices
                 'String','Sens. Cutoff',...
                 'Position',[sensx+10,sensy,100,20],... % position on UI
                 'FontName','Arial','FontWeight','bold','FontSize',12,...
                 'BackgroundColor',[0.8,0.8,0.8]};
UI(2).Control = {'Style','slider',...        % popupmenu gives list of choices
                 'Position',[sensx,sensy-30,120,30],... % position on UI
                 'Max',1.0,'Min',0,'Value',Recon(1).senscutoff,...
                 'SliderStep',[0.025 0.1],...
                 'Callback',{@sensCutoffCallback}};
UI(2).Callback = {'sensCutoffCallback.m',...
                 'function sensCutoffCallback(hObject,eventdata)',...
                 ' ',...
                 'sens = get(hObject,''Value'');',...
                 'ReconL = evalin(''base'', ''Recon'');',...
                 'for i = 1:size(ReconL,2)',...
                 '    ReconL(i).senscutoff = sens;',...
                 'end',...
                 'assignin(''base'',''Recon'',ReconL);',...
                 '% Set Control.Command to re-initialize Recon structure.',...
                 'Control = evalin(''base'',''Control'');',...
                 'Control.Command = ''update&Run'';',...
                 'Control.Parameters = {''Recon''};',...
                 'assignin(''base'',''Control'', Control);',...
                 '% Set the new cutoff value in the text field.',...
                 'h = findobj(''tag'',''sensCutoffValue'');',...
                 'set(h,''String'',num2str(sens,''%1.3f''));',...
                 'return'};
UI(3).Control = {'Style','edit','String',num2str(Recon(1).senscutoff,'%1.3f'), ...  % text
                 'Position',[sensx+20,sensy-40,60,22], ...   % position on UI
                 'tag','sensCutoffValue', ...
                 'BackgroundColor',[0.9,0.9,0.9]}; 
     
 % -- Enable DisplayWindow's WindowButtonDown callback function for switching acquisition loops.
UI(4).Control = {'UserC2', 'Style', 'VsPushButton', 'Label', 'Trigger Acq.'};
UI(4).Callback = text2cell('%TriggeredAcquisition');

clear i j n sensx sensy

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 2;

% Save all the structures to a .mat file.
save(['./MatFiles/' scriptName '.mat']);
display(['filename =''' scriptName '.mat'';VSX'])
% eval(['filename =''' scriptName ''';VSX'])


return

%TriggeredAcquisition
Control =evalin('base', 'Control');
Control(1).Command='set';
lastBmodeEvent = evalin('base', 'lastBmodeEvent');
Control(1).Parameters = {'Parameters', 1, 'startEvent', lastBmodeEvent};
evalin('base', sprintf('Resource.Parameters.startEvent = %d;', lastBmodeEvent));
assignin('base','Control', Control);
disp('Switching to triggered acquisition')
return
%TriggeredAcquisition
