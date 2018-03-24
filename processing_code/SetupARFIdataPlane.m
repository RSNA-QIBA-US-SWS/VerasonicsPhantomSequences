
%==========================================================================

function [dispPlane,tms,latmm] = SetupARFIdataPlane(arfidata,axialmm,latmm,tms,pushFocusMM,fliplat,par)

        %    

    DOFmm           = par.DOFmm;
    minLatMM        = par.minLatMM;
    maxLatMM        = par.maxLatMM;
    maxTimeMS       = par.maxTimeMS;
    minTimeMS       = par.minTimeMS;
    LPFcutoffKHz    = par.LPFcutoffKHz;
    nStepsToRemove  = par.nStepsToRemove;
    desiredPRFkHz   = par.desirecPRFkHz;
    rolloffTimeMS   = par.rolloffTimeMS;

    axialIdx = find(axialmm>=pushFocusMM-DOFmm/2 & axialmm<=pushFocusMM+DOFmm/2);
    dispPlane = squeeze(mean(arfidata(axialIdx,:,:)));

    if fliplat
        dispPlane = flipdim(dispPlane,1);   % no flip for latmm since symmetric
    end

    latidx = find(latmm>=minLatMM & latmm <= maxLatMM);     % limit lateral range
    dispPlane = dispPlane(latidx,:);
    latmm = latmm(latidx);

    [dispPlane,tms,latmm] = RemoveAndFixReverb(dispPlane,tms,latmm,nStepsToRemove);

    dispPlane = LowPassFilterInTime2D(dispPlane,tms,LPFcutoffKHz);

    tIdx = find(tms<=maxTimeMS);        % truncate to maxTime
    dispPlane = dispPlane(:,tIdx);
    tms = tms(tIdx);

    [dispPlane,tms] = AddZerosBeforePush(dispPlane,tms,minTimeMS);

    [dispPlane,tms] = UpsampleTimeAndDispPlane(dispPlane,tms,desiredPRFkHz);

    dispPlane=RollOffPlaneTimeEdges(dispPlane,tms,rolloffTimeMS);
end

%==========================================================================

function [plane,tms,latmm] = RemoveAndFixReverb(plane,tms,latmm,nStepsToRemove)

    tidx = find(sum(abs(plane))==0);    % last reference has disp=0
    if length(tidx)~=1
        error('error')
    end
    plane = plane(:,tidx+1:end);        % remove timesteps before push
    tms = tms(tidx+1:end);

    plane = plane(:,nStepsToRemove+1:end);
    tms = tms(nStepsToRemove+1:end);

    for ilat=1:length(latmm)
        plane(ilat,:) = plane(ilat,:) - plane(ilat,1);
    end
end

%==========================================================================

function output = LowPassFilterInTime2D(plane,t,fcutoff)

    [nlats,ntimes]=size(plane);    
    output=zeros(nlats,ntimes);

    fs=1/mean(diff(t));         % sampling frequency
    fnorm=fcutoff/(fs/2);
    
    filtorder=3;
    [b,a]=butter(filtorder,fnorm,'low');

    for ilat=1:nlats
        output(ilat,:)=filtfilt(double(b),double(a),plane(ilat,:));
    end
end

%==========================================================================

function [plane,tms] = AddZerosBeforePush(plane,tms,minTimeMS)

    dt=mean(diff(tms));                     % determine number of steps to add
    nsteps = floor((tms(1)-minTimeMS)/dt);

    tToAdd = tms(1) - (nsteps:-1:1)*dt;         % determine new timesteps
    planeToAdd = zeros(size(plane,1),nsteps);   %   and make plane of zeros

    tms = [tToAdd(:);tms(:)];
    plane = [planeToAdd plane];
end

%==========================================================================

function [planeup,tup]=UpsampleTimeAndDispPlane(plane,t,desiredPRF)

    expPRF = 1/mean(diff(t));
    upsamp = round(desiredPRF/expPRF);

    deltat = mean(diff(t))/upsamp;
    tup=t(1):deltat:max(t);

    [nlats,ntimes]=size(plane);
    planeup=zeros(nlats,length(tup));

    for ilat=1:nlats
        planeup(ilat,:)=interp1(t,plane(ilat,:),tup,'spline');
    end
end

%==========================================================================

function newplane=RollOffPlaneTimeEdges(plane,tms,rolloffTimeMS)

    [nlats,ntimes]=size(plane);
    tmask=ones(1,ntimes);
    twidth = rolloffTimeMS/mean(diff(tms));

    for itime=1:ntimes
        if itime < (twidth+1)
            tmask(itime) = 0.5*(1-cos( (itime-1)/twidth*pi ));
        end
        if itime > ntimes-twidth
            tmask(itime) = 0.5*(1+cos((itime-(ntimes-twidth))/twidth*pi));
        end
    end

    mask=ones(size(plane));
    for ilat=1:nlats
        mask(ilat,:)=mask(ilat,:).*tmask;
    end
    newplane = plane.*mask;
end

%==========================================================================
