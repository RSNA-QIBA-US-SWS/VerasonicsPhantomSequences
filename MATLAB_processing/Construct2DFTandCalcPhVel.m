
%==========================================================================

function phVel = Construct2DFTandCalcPhVel(plane,t,lat,freqsToAnalyzeHz)

        %

    nt=2000;                % t in sec, ft in Hz
    dt=mean(diff(t));
    fs=1/dt;
    ft=(-nt/2:nt/2-1)*(fs/nt);

    nx=2000;                % lat in m, flat in cyc/m
    dlat=mean(diff(lat));
    fslat=1/dlat;
    flat=(-nx/2:nx/2-1)*(fslat/nx);

    absfftplane = abs(fftshift(fft2(plane,nx,nt)));

    phVel = NaN(1,length(freqsToAnalyzeHz));
    for ifreq=1:length(freqsToAnalyzeHz)

        ftidx = find(roundn(ft,-2) == roundn(freqsToAnalyzeHz(ifreq),-2));
        if length(ftidx)~=1
            error('error finding ftidx')
        end        
        profile = absfftplane(:,ftidx);
        [maxval flatidx] = max(profile);   % find max along fixed freq
        
        if flatidx==1 | flatidx==length(profile)
            maxspfval = flat(flatidx);        % maxk is spatial frequency of maximum signal
        else
            y1 = profile(flatidx-1);       % interpolate for better estimate
            y2 = profile(flatidx);
            y3 = profile(flatidx+1);
            interpval = (y1-y3)/(y1-2*y2+y3)/2;
            maxspfval = flat(flatidx) + mean(diff(flat))*interpval;
        end

        phVel(ifreq) = -freqsToAnalyzeHz(ifreq)/maxspfval;        % phase velocity from peak location
    end
end

%==========================================================================
