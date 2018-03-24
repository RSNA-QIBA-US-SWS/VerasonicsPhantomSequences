function genDispMTL(filestamp)

if nargin < 1
paramFile = dir('*_parameters.mat');
filestamp = paramFile(end).name(1:14);
end

if isnumeric(filestamp)
    paramFile = dir('*_parameters.mat');
    filestamp = paramFile(filestamp).name(1:14);
end

paramfile = [filestamp '_parameters.mat'];
S = load(paramfile,'Vvalue','PData','TX','Resource','Trans','ne','nrefs','lastBmodeTransmit','T','T_idx','TW');
Vvalue              = S.Vvalue;
PData               = S.PData;
TX                  = S.TX;
TW                  = S.TW;
Resource            = S.Resource;
Trans               = S.Trans;
lastBmodeTransmit   = S. lastBmodeTransmit;
ne                  = S.ne;
nrefs               = S.nrefs;
T                   = S.T;
T_idx               = S.T_idx;
clear S

c = Resource.Parameters.speedOfSound;
fc= Trans.frequency*1e6;
fTrack = TW(1).Parameters(1);
if fTrack>10
    fTrack = (180/fTrack/2).*1e6;
else
    fTrack = fTrack.*1e6;
end
push_focus = TX(lastBmodeTransmit+2).focus*c/fc;
% V = [Vvalue(1) Vvalue(end)];

lat = ([0:PData(2).Size(2)-1].*PData(2).PDelta(1)+PData(2).Origin(1)).*c/fc*100;
axial = ([0:PData(2).Size(1)-1].*PData(2).PDelta(3)+PData(2).Origin(3)).*c/fc*100;

idx = 1;
nvalstoread = PData(2).Size(1)*PData(2).Size(2)*ne;
ptr = nvalstoread*(idx-1)*4;
fid = fopen([filestamp '_IQreal.bin'],'rb');
fseek(fid,ptr,'bof');
I = fread(fid,nvalstoread,'int32');
fclose(fid);

fid = fopen([filestamp '_IQimag.bin'],'rb');
fseek(fid,ptr,'bof');
Q = fread(fid,nvalstoread,'int32');
fclose(fid);

IQ = I + 1j*Q;
IQ = reshape(IQ,PData(2).Size(1),PData(2).Size(2),ne);
dimIQ = size(IQ);

k_length = floor(1.5*4/PData(2).PDelta(3)/2)*2+1; % 1.5 wavelength
arfidata = -computeDisplacementsSamTrack(IQ,nrefs-1,k_length,c,fTrack);
[arfidata]=phase_unwrap(arfidata*1e6,fTrack,nrefs+2,1,0);
focalPush = round(push_focus(ceil(idx/2.5))*1000);
% voltage = V(2-mod(idx,2));

% keyboard

% [dispout] = motion_filter(arfidata,T,T_idx,2);

numrefs = nrefs;
arfidata = single(arfidata);
axial = single(axial);
lat = single(lat);
save([filestamp '_fromIQ_arfidata'],'arfidata','axial','lat','numrefs','push_focus','T')


end
