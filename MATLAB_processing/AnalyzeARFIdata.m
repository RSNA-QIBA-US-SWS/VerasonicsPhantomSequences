
%==========================================================================

function AnalyzeARFIdata(filestamp,par)

    %

    filename = [par.analysisDir '/' filestamp '_fromIQ_arfidata.mat'];
    S = load(filename,'arfidata','axial','lat','numrefs','push_focus','T');
    arfidata      = double(S.arfidata);           % arfidata saved in um
    axialmm       = double(S.axial) *10;          % axial saved in cm
    latmmARFIdata = double(S.lat) *10;            % lat saved in cm
    pushFocusMM   = double(S.push_focus) *1000-2;   % push_focus saved in m
    tmsARFIdata   = double(S.T);                  % T saved in ms
    
    clear S

    for fliplat=0:1

        [dispPlane,tms,latmm] = SetupARFIdataPlane(arfidata,axialmm,latmmARFIdata,tmsARFIdata,pushFocusMM,fliplat,par);

        velPlane = diff(dispPlane,[],2)./mean(diff(tms)); % diff. for particle vel plane.
        veltms = (tms(1:end-1)+tms(2:end))/2;
        vellatmm = latmm;

        TTPdisp = [];TTPvel = [];
        
        if ismember(lower(par.swsest),{'ttp',''})
            [TTPdisp,Vdisp] = FindTTPandGSWS(dispPlane,tms,   latmm);       % get gSWSs using TTP
            [TTPvel, Vvel]  = FindTTPandGSWS(velPlane, veltms,vellatmm);
        end

        wtVelPlane = NaN(size(velPlane)); % apply sqrt(x) weighting to velPlane
        for it=1:length(veltms)
            wtVelPlane(:,it) = sqrt(vellatmm(:)).*velPlane(:,it);
        end

        t = veltms /1000;  % convert veltms, vellatmm to MKS
        lat = vellatmm / 1000;
        freqsToAnalyzeHz = par.freqsToAnalyzeHz;
        phVel = Construct2DFTandCalcPhVel(wtVelPlane,t,lat,freqsToAnalyzeHz);

        % plotting function
        if par.plotfig
            figure(99);
            subplot(1,2,1)
            imagesc(tms,latmm,dispPlane);hold on;
            scatter(TTPdisp, latmm, 'filled', 'k');
            xlabel('Time (ms)');ylabel('Lateral (mm)');title('Particle Displacement');
            ylim([min(latmm) max(latmm)]);
    
            subplot(1,2,2)
            imagesc(tms,latmm,velPlane);hold on;
            scatter(TTPvel, latmm, 'filled', 'k');
            xlabel('Time (ms)');ylabel('Lateral (mm)');title('Particle velocity');
            ylim([min(latmm) max(latmm)]);
            pause(0.5)
        end

        saveFile = [par.analysisDir '/' filestamp '_FL' num2str(fliplat) '_phVel_gSWS_data.mat'];
        save(saveFile,'dispPlane','tms','latmm','velPlane','veltms','vellatmm','TTPdisp','TTPvel','Vdisp','Vvel','freqsToAnalyzeHz','phVel')
    end
end

%==========================================================================

function [ttp,gSWS] = FindTTPandGSWS(plane,tms,latmm)

    ttp = NaN(1,length(latmm));
    for ilat=1:length(latmm)
        [maxval, maxidx] = max(plane(ilat,:));
        ttp(ilat) = tms(maxidx);
    end

    p = polyfit(latmm,ttp,1);
    gSWS = 1/p(1);
end

% =========================================================================

function [latsums]=make_latsum(plane,itminstart,itmaxstart,itmaxend)

[nlats,ntimes]=size(plane);

interpIndx=zeros(ntimes-1,nlats);         % get indices and fractions for
interpFrac=zeros(ntimes-1,nlats);         % interpolation at lat positions

ilats=1:nlats;

for idiff=1:ntimes-1   %run across every time step
    tvals=(idiff-1)/(length(ilats)-1)*(ilats-ilats(1));  % 'ideal' time steps across each lat to connect the way we'd want to
    interpIndx(idiff,ilats)=floor(tvals);
    interpFrac(idiff,ilats)=1-(tvals-interpIndx(idiff,ilats));
end

latsums = nan(ntimes-1,ntimes-1);
maxistart=min(itmaxstart,ntimes-1); %use given index, or last index available if smaller
maxiend=min(itmaxend,ntimes-1);
startTimeIndex=itminstart;
endoffset = 0;

for istart=startTimeIndex:maxistart
    iendvector=istart+endoffset:maxiend;  
    idiff=iendvector-istart+1;            % +1 for matlab numbering  Find difference between beginning and end time steps
    
    idx0vals=interpIndx(idiff,ilats)+istart;  %Index zero vals
    fracs=interpFrac(idiff,ilats);
    
    indices1=ones(size(iendvector))'*ilats + (idx0vals-1)*size(plane,1); %Each 1x28 matrices of indexes-- drawing the 'line'
    indices2=ones(size(iendvector))'*ilats+ (idx0vals)*size(plane,1); % (id0vals+1)-1
    latsums(istart,iendvector) = sum(plane(indices1).*fracs+plane(indices2).*(1-fracs),2,'omitnan')';  % sum over lat pos,weighting the data points on either side of the 'line' appropriately
end
end
