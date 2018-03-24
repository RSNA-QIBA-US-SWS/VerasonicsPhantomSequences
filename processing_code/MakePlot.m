
%==========================================================================

function MakePlot

        % plot parameters

    titleStr     = 'E2297-A3';      % title string
    maxPlotSpeed = 4.5;             % maximum speed for plots
    maxplotfreq  = 820;             % max freq to plot
    fracErrorBar = 0.3;             % fixed fraction for error bars
    CIfactor     = 1.96;            % confidence interval for error bars
    fmax         = 800;             % maximum freq to plot
    rmin         = 0.004;           % minimum lateral position for low freq cutoff
    krThreshold  = 1.5;             % kr threshold for low freq cutoff
    maxSpeed     = 8;               % max speed for "good" result
    thFactor     = 2;               % phVel accept range = +/- thfactor * std

    freqsToAnalyzeHz = 20:10:800;   % plot frequencies

    paramFiles = dir('*_parameters.mat');       % get parameter files

    Vdisp = NaN(1,2*length(paramFiles));
    Vvel  = NaN(1,2*length(paramFiles));
    phVel = NaN(2*length(paramFiles),length(freqsToAnalyzeHz));

    for ifile = 1:length(paramFiles)                % load Vdisp, Vvel, phVel data
        filestamp = paramFiles(ifile).name(1:14);

        for fliplat=0:1
            idx = ifile*2+fliplat-1;
            dataFile = [filestamp '_FL' num2str(fliplat) '_phVel_gSWS_data.mat'];
            S=load(dataFile,'Vdisp','Vvel','freqsToAnalyzeHz','phVel');
            Vdisp(idx)       = S.Vdisp;
            Vvel(idx)        = S.Vvel;
            phVel(idx,:)     = S.phVel;
            freqsToAnalyzeHz = S.freqsToAnalyzeHz;
            clear S
        end
    end

    Vdisp_avg = nanmean(Vdisp);     % find avg,std gSWS results
    Vdisp_std =  nanstd(Vdisp);
    Vvel_avg  = nanmean(Vvel);
    Vvel_std  =  nanstd(Vvel);

    [phVel_avg,phVel_std] = MeanStdOmitBadPtsOneCase_usingPositiveData(phVel,thFactor,maxSpeed);

    black = [0.0 0.0 0.0];
    gray  = [0.4 0.4 0.4];
    x1=1.0;
    x2=2.0;
    dx = 0.04;
    YTickVals = [0:1:maxPlotSpeed];

    figure(1)
    clf
    set(gcf,'Units','inches','Position',[0.5 2 8 3.5])

    hhh1=subplot(1,2,1);
    hold on
    ybot = (1-fracErrorBar)*Vdisp_avg;      % +/- fraction of avg displacement gSWS
    ytop = (1+fracErrorBar)*Vdisp_avg;
    plot([x1 x1], [ybot ytop],      'Color',gray, 'LineWidth',1)
    plot([x1-dx x1+dx],[ybot ybot], 'Color',gray, 'LineWidth',1)
    plot([x1-dx x1+dx],[ytop ytop], 'Color',gray, 'LineWidth',1)

    ybot = Vdisp_avg - CIfactor*Vdisp_std;          % +/- CI
    ytop = Vdisp_avg + CIfactor*Vdisp_std;
    plot([x1 x1], [ybot ytop],      'Color',black,'LineWidth',2)
    plot([x1-dx x1+dx],[ybot ybot], 'Color',black,'LineWidth',2)
    plot([x1-dx x1+dx],[ytop ytop], 'Color',black,'LineWidth',2)

    ybot = (1-fracErrorBar)*Vvel_avg;       % +/- 30% of avg velocity gSWS
    ytop = (1+fracErrorBar)*Vvel_avg;
    plot([x2 x2], [ybot ytop],      'Color',gray, 'LineWidth',1)
    plot([x2-dx x2+dx],[ybot ybot], 'Color',gray, 'LineWidth',1)
    plot([x2-dx x2+dx],[ytop ytop], 'Color',gray, 'LineWidth',1)

    ybot = Vvel_avg - CIfactor*Vvel_std;            % +/- CI
    ytop = Vvel_avg + CIfactor*Vvel_std;
    plot([x2 x2], [ybot ytop],      'Color',black,'LineWidth',2)
    plot([x2-dx x2+dx],[ybot ybot], 'Color',black,'LineWidth',2)
    plot([x2-dx x2+dx],[ytop ytop], 'Color',black,'LineWidth',2)
    hold off

    set(gca,'Units','inches','YTick',YTickVals)
    set(gca,'XTick',[1 2],'XTickLabel',{'disp' 'vel'})
    set(gca,'Box','on')
    ylim([0 maxPlotSpeed])
    xlim([0 3])
    ylabel('group SWS (m/s)')
    title(titleStr)
    set(gca,'Units','inches','Position',[0.5 0.5 1.4 2.7])

    kr = 2*pi*freqsToAnalyzeHz./phVel_avg.*rmin;
    fidx = find(kr>=krThreshold & freqsToAnalyzeHz<=fmax);
    ddf=2;
    df=10;
    nn=length(fidx);
    patchx = NaN(1,2*nn+2);
    patchy = NaN(1,2*nn+2);

    for ii=1:nn+1
        switch ii
            case 1
                patchx(ii)            = freqsToAnalyzeHz(fidx(ii))-df/2;
                patchx(2*(nn+1)-ii+1) = freqsToAnalyzeHz(fidx(ii))-df/2;
                patchy(ii)            = (1+fracErrorBar)*phVel_avg(fidx(ii));
                patchy(2*(nn+1)-ii+1) = (1-fracErrorBar)*phVel_avg(fidx(ii));
            case nn+1
                patchx(ii)            = freqsToAnalyzeHz(fidx(ii-1))+df/2;
                patchx(nn+2)          = freqsToAnalyzeHz(fidx(ii-1))+df/2;
                patchy(ii)            = (1+fracErrorBar)*phVel_avg(fidx(ii-1));
                patchy(nn+2)          = (1-fracErrorBar)*phVel_avg(fidx(ii-1));
            otherwise
                patchx(ii)            = freqsToAnalyzeHz(fidx(ii))-df/2;
                patchx(2*(nn+1)-ii+1) = freqsToAnalyzeHz(fidx(ii))-df/2;
                patchy(ii)            = (1+fracErrorBar)/2*(phVel_avg(fidx(ii-1))+phVel_avg(fidx(ii)));
                patchy(2*(nn+1)-ii+1) = (1-fracErrorBar)/2*(phVel_avg(fidx(ii-1))+phVel_avg(fidx(ii)));
        end
    end    

    hhh2=subplot(1,2,2);
    hold on
    fill(patchx,patchy,[0.8 0.8 0.8],'EdgeColor','none');

    for ii=1:length(fidx)

        ybot = phVel_avg(fidx(ii)) - CIfactor*phVel_std(fidx(ii));        % +/- CI
        ytop = phVel_avg(fidx(ii)) + CIfactor*phVel_std(fidx(ii));
        plot([freqsToAnalyzeHz(fidx(ii))     freqsToAnalyzeHz(fidx(ii))    ],[ybot ytop], 'Color',black,'LineWidth',2)
        plot([freqsToAnalyzeHz(fidx(ii))-ddf freqsToAnalyzeHz(fidx(ii))+ddf],[ybot ybot], 'Color',black,'LineWidth',2)
        plot([freqsToAnalyzeHz(fidx(ii))-ddf freqsToAnalyzeHz(fidx(ii))+ddf],[ytop ytop], 'Color',black,'LineWidth',2)
    end
    hold off
    set(gca,'Units','inches','YTick',YTickVals)
    set(gca,'Box','on')
    ylim([0 maxPlotSpeed])
    xlim([0 maxplotfreq])
    ylabel('phase velocity (m/s)')
    xlabel('frequency (Hz)')
    title(titleStr)
    set(gca,'Units','inches','Position',[2.4 0.5 5.4 2.7])

    savefig('gSWS_phVel_figure')
end

%==========================================================================

function [out_avg,out_std] = MeanStdOmitBadPtsOneCase_usingPositiveData(data,factor,maxSpeed)

    [npts,nfreqs] = size(data);
    out_avg = NaN(1,nfreqs);
    out_std = NaN(1,nfreqs);

    for ifreq=1:nfreqs
        temp = data(:,ifreq);
        idx = find(temp>0 & temp<maxSpeed);
        temp_avg = mean(temp(idx));
        temp_std =  std(temp(idx));

        out_avg(ifreq) = temp_avg;      % set outputs in case something goes wrong
        out_std(ifreq) = temp_std;      %   with goodidx below
        
        if temp_std>0
            goodidx = find( abs(temp-temp_avg)<factor*temp_std & temp>0);        
        else
            goodidx=1:npts;
        end
        out_avg(ifreq) = nanmean(temp(goodidx));
        out_std(ifreq) =  nanstd(temp(goodidx));
    end
end

%==========================================================================
