
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
        if par.plot_int_fig
            figure(99);
            subplot(1,2,1)
            imagesc(tms,latmm,dispPlane);hold on;
            scatter(TTPdisp, latmm, 'filled', 'k');
            xlabel('Time (ms)');ylabel('Lateral (mm)');title('Particle Displacement');
            ylim([min(latmm) max(latmm)]);
            axis image
    
            subplot(1,2,2)
            imagesc(tms,latmm,velPlane);hold on;
            scatter(TTPvel, latmm, 'filled', 'k');
            xlabel('Time (ms)');ylabel('Lateral (mm)');title('Particle velocity');
            ylim([min(latmm) max(latmm)]);
            axis image
            pause(0.2)
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