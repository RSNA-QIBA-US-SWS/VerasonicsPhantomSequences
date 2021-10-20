clear all
global filedir
scriptName='SETUPC5_2SHEAR_WAVE_MTL';

% USER MUST EDIT THIS PATH
filedir = '/home/verasonics/cloud/Vantage-4.2.0-2001220500';
cd(filedir)
push_focus = 50; %mm
push_Fnum = 1.5; % 
pushCycle  = 900;
ne = 205;
nrefs = 5;
pushAngleDegree = 0;

%% Define basic parameters
m = 64; % Bmode lines
bmode_focus = 50;
bmode_Fnum = 2;

getPower = 0; % DO NOT delete, used by save_swei_data %----------------------%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
saveChannelData = 0;
saveIQData = 1;

Resource.Parameters.connector = 1;
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.simulateMode = 0;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;

%% Specify Trans structure array.
Trans.name = 'C5-2';
Trans.units = 'wavelengths';
Trans = computeTrans(Trans);  % C5-2 transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 76;  % set maximum high voltage limit for pulser supply.
TPC(5).maxHighVoltage = 76;
radius = Trans.radius;
w = Resource.Parameters.speedOfSound/Trans.frequency/1000; % wavelength in mm

%% Specify angles
scanangle = Trans.numelements*Trans.spacing/radius;
aperture_size = radius*scanangle; % mm
Angle = (-scanangle/2):(scanangle/(m-1)):(scanangle/2); % Bmode angle range
dtheta = scanangle/(m-1);

P.numRays = m;  % no. of raylines to program
P.startDepth = 0; % startDepth and endDepth in wavelength
P.endDepth = 208;

%% Specify PData structure array.
PData(1).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5; 
PData(1).Origin(1,2) = 0;
PData(1).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;
% Define PData Regions for numRays scanlines
for j = 1:P.numRays
    PData(1).Region(j) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[0,0,-radius],...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta,...
                       'steer',Angle(j)));
end
PData(1).Region = computeRegions(PData(1));


% Specify PData structure array.
PData(2).PDelta = [0.2/w, 0, 0.5];  % x, y and z pdeltas
sizeRows2 = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/4)))/PData(2).PDelta(3));
sizeCols2 = 9 + ceil(2*(P.endDepth + radius)*sin(scanangle/4)/PData(2).PDelta(1));
PData(2).Size = [sizeRows2,sizeCols2,1];     % size, origin and pdelta set region of interest.
PData(2).Origin(1,1) = -floor(sizeCols2/2)*0.2/w; 
PData(2).Origin(1,2) = 0;
PData(2).Origin(1,3) = ceil(radius * cos(-scanangle/4)) - radius - 5;

PData(2).Region = struct(...
    'Shape',struct('Name','Sector','Position',[0,0,-radius],'r1',radius+P.startDepth,'r2',radius+P.endDepth,'angle',scanangle/2));


%% Specify resource buffers
% RcvBuffer stores channel data. Buffer 1 stores Bmode data
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = 2400*m; % max 4096 per axial line
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 2;  

Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = ne*2500;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 1;

Resource.VDAS.dmaTimeout = 8000;

Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
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
% - Bmode and track waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.8,2,1];   % A, B, C, D

% - Push waveform.
TW(2).type = 'parametric';
TW(2).Parameters = [2.3585,1,pushCycle*2,1];  %

%% Specify TX structure array.  
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', round(bmode_focus/w), ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, m+2);
               
% B-Mode TX's               
onele_bmode = round(((bmode_focus/bmode_Fnum)/(aperture_size*w)*Resource.Parameters.numRcvChannels)/4)*4;
offele_bmode = Resource.Parameters.numRcvChannels-onele_bmode;
for n = 1:m
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    TX(n).Apod = [zeros(1,ceil(n/m*offele_bmode)) ones(1,onele_bmode) zeros(1,offele_bmode-ceil(n/m*offele_bmode))];
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
n = m+2;
onele_push = round(((push_focus/push_Fnum)/(aperture_size*w)*Resource.Parameters.numRcvChannels)/4)*4;
offele_push = Resource.Parameters.numRcvChannels-onele_push;
TX(n).waveform = 2;
TX(n).Origin = [0.0 0.0 0.0];
TX(n).focus = round(push_focus/w);
TX(n).Apod = [zeros(1,offele_push/2) ones(1,onele_push) zeros(1,offele_push/2)];
TX(n).Delay = computeTXDelays(TX(n));

%% Specify Receive structure arrays. 
% Compute the maximum receive path length, using the law of cosines.
maxAcqLength = sqrt((P.endDepth+radius)^2 + radius^2 - ...
                     2*(P.endDepth+radius)*radius*cos(scanangle)) - P.startDepth;
wlsPer128 = 128/(4*2); % wavelengths in 128 samples for 4 samplesPerWave
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wlsPer128*ceil(maxAcqLength/wlsPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, Resource.RcvBuffer(1).numFrames*m+ne);

for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(m*(i-1)+1).callMediaFunc = 1;
    for j = 1:m
        Receive(m*(i-1)+j).framenum = i;
        jj=2*j;
        lft = jj - 32;
        if lft < 1, lft = 1; end;
        rt = jj + 32;
        if rt > Trans.numelements, rt = Trans.numelements; end;
        Receive(m*(i-1)+j).Apod(lft:rt) = 1.0;
        Receive(m*(i-1)+j).acqNum = j; 
    end
end

lastBmodeReceive = m*(i-1)+j;  

% - Set event specific Receive attributes for SWEI
for j = 1:ne
    Receive(lastBmodeReceive+j).Apod(:) = 1.0;
    Receive(lastBmodeReceive+j).bufnum = 2;
    Receive(lastBmodeReceive+j).framenum = 1;
    Receive(lastBmodeReceive+j).acqNum = j;
    Receive(lastBmodeReceive+j).startDepth = 0;
end

%% Specify TGC Waveform structure.
TGC.CntrlPts = [500,590,650,710,770,800,850,950];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

%% Specify Recon structure arrays.
Recon(1) = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'newFrameTimeout',1000,...
               'IntBufDest', [0,0], ...
               'ImgBufDest', [1,-1], ...
               'RINums', [1:m]);
%                'rcvBufFrame',[1,-1], ...

Recon(2) = struct('senscutoff', 0.6, ...
               'pdatanum', 2, ...
               'newFrameTimeout',2000,...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [0,0], ...
               'RINums',(m+1:(m+ne))');

%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 0, ...  % replace data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 0), 1, m);
% - Set specific ReconInfo attributes.
for i = 1:m
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end

% - ReconInfo for SWEI.
k = m; 
ReconInfo((k+1):(k+ne)) = repmat(struct('mode', 3, ... % IQ output 
                   'txnum', m+1, ...
                   'rcvnum', lastBmodeReceive+1, ...
                   'regionnum', 1), 1, ne);

% - Set specific ReconInfo attributes.
for j = 1:ne 
    ReconInfo(k+j).rcvnum = lastBmodeReceive+j;
    ReconInfo(k+j).pagenum = j;
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
                         'displayWindow',1};
                            
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
                         'srcframenum',0,... %1  -- 0 is most recent, which will be 1 in this case, as there is only 1 frame in the buffer
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
SeqControl(10).argument = 400;
% - time between last Bmode tx and jump back 1st Bmode tx to include time
SeqControl(11).command = 'timeToNextAcq';
SeqControl(11).argument = 500500;
% - wait for transfer before marking processed
% - Jump back to start.
SeqControl(12).command = 'jump';
SeqControl(12).argument = 1;
% - time between 'acquisitions', not used

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

lastBmodeEvent = n

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

% SWEI loop
%push 1 ensemble

for j = 1:nrefs
    Event(n).info = 'Acquire reference data';
    Event(n).tx = m+1;
    Event(n).rcv = lastBmodeReceive+j;
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
Event(n).tx = m+2;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 10;
TTNAS(T_i)=TTNAS(T_i-1)+SeqControl(Event(n).seqControl).argument;
n = n+1;

for j = nrefs+1:ne
    Event(n).info = 'Acquire tracking data';
    Event(n).tx = m+1;
    Event(n).rcv = lastBmodeReceive+j;
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

%         Event(n).info = 'Wait for transfer complete';
%         Event(n).tx = 0;         % no transmit
%         Event(n).rcv = 0;        % no rcv
%         Event(n).recon = 0;      % reconstruction
%         Event(n).process = 0;    % process
%         Event(n).seqControl = nsc; % waitForTransferComplete
%         SeqControl(nsc).command = 'waitForTransferComplete'; % transfer frame to host buffer
%         SeqControl(nsc).argument = nsc-1;
%         nsc = nsc+1;
%         n = n+1;

Event(n).info = 'recon';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 2;      % reconstruction
Event(n).process = 0;    % process
Event(n).seqControl = 2; % NOOP for .1sec
n = n+1;

Event(n).info = 'noop for charging ext. cap.';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 2;
n = n+1;

% if saveChannelData
%     Event(n).info = 'save channel data';
%     Event(n).tx = 0;
%     Event(n).rcv = 0;
%     Event(n).recon = 0;
%     Event(n).process = 2;
%     Event(n).seqControl = 2;
%     n = n+1;
% end

Event(n).info = 'save IQ data';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 3;
Event(n).seqControl = 2;
n = n+1;

Event(n).info = 'Back to Matlab';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % no process
Event(n).seqControl = 5; % Back to Matlab
n = n+1;

Event(n).info = 'Jump back to third event (stay at current power)';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
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
