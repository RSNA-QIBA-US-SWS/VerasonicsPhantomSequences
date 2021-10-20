
%==========================================================================

function errorFlag = MakePlotAndSaveOutputCSVfile(par)

        %

    analysisDir      = par.analysisDir;         % directory with result files saved in analysis
    saveDir          = par.saveDir;             % directory to save output plot and CSV data
    
    freqsToAnalyzeHz = par.freqsToAnalyzeHz;    % analysis frequencies (Hz)
    phantomID        = par.phantomID;
    maxPlotSpeed     = par.maxPlotSpeed;
    maxplotfreq      = par.maxplotfreq;
    fracErrorBar     = par.fracErrorBar;
    CIfactor         = par.CIfactor;
    fmax             = par.fmax;
    rmin             = par.rmin;
    krThreshold      = par.krThreshold;
    maxSpeed         = par.maxSpeed;
    thFactor         = par.thFactor;


    analysisFiles = dir([analysisDir '/*_FL*_phVel_gSWS_data.mat']);    % get analysis files

    Vdisp = NaN(1, length(analysisFiles));
    Vvel  = NaN(1, length(analysisFiles));
    phVel = NaN(2 * length(analysisFiles), length(freqsToAnalyzeHz));

    for ifile = 1:length(analysisFiles)         % load Vdisp, Vvel, phVel data

        S = load(analysisFiles(ifile).name, 'Vdisp', 'Vvel', 'freqsToAnalyzeHz', 'phVel');
        Vdisp(ifile)     = S.Vdisp;
        Vvel(ifile)      = S.Vvel;
        phVel(ifile, :)  = S.phVel;
        freqsToAnalyzeHz = S.freqsToAnalyzeHz;
        clear S
    end

    [Vdisp_avg, Vdisp_std] = MeanStdOmitBadPtsOneCase_gSWS_usingPositiveData(Vdisp, thFactor, maxSpeed);
    [Vvel_avg, Vvel_std]   = MeanStdOmitBadPtsOneCase_gSWS_usingPositiveData(Vvel, thFactor, maxSpeed);

    [phVel_avg, phVel_std] = MeanStdOmitBadPtsOneCase_usingPositiveData(phVel, thFactor, maxSpeed);

    black = [0.0 0.0 0.0];      % make plot
    gray  = [0.4 0.4 0.4];
    x1 = 1.0;
    x2 = 2.0;
    dx = 0.04;
    YTickVals = 0:1:maxPlotSpeed;

    figure(1)
    clf
    set(gcf, 'Units', 'inches', 'Position', [0.5 2 8 3.5])

    hhh1 = subplot(1, 2, 1);
    hold on
    ybot = (1 - fracErrorBar) * Vdisp_avg;      % +/- fraction of avg displacement gSWS
    ytop = (1 + fracErrorBar) * Vdisp_avg;
    plot([x1 x1], [ybot ytop],       'Color',gray, 'LineWidth',1)
    plot([x1-dx x1+dx], [ybot ybot], 'Color',gray, 'LineWidth',1)
    plot([x1-dx x1+dx], [ytop ytop], 'Color',gray, 'LineWidth',1)

    ybot = Vdisp_avg - CIfactor * Vdisp_std;          % +/- CI
    ytop = Vdisp_avg + CIfactor * Vdisp_std;
    plot([x1 x1], [ybot ytop],       'Color',black,'LineWidth',2)
    plot([x1-dx x1+dx], [ybot ybot], 'Color',black,'LineWidth',2)
    plot([x1-dx x1+dx], [ytop ytop], 'Color',black,'LineWidth',2)

    ybot = (1 - fracErrorBar) * Vvel_avg;       % +/- 30% of avg velocity gSWS
    ytop = (1 + fracErrorBar) * Vvel_avg;
    plot([x2 x2], [ybot ytop],       'Color',gray, 'LineWidth',1)
    plot([x2-dx x2+dx], [ybot ybot], 'Color',gray, 'LineWidth',1)
    plot([x2-dx x2+dx], [ytop ytop], 'Color',gray, 'LineWidth',1)

    ybot = Vvel_avg - CIfactor * Vvel_std;            % +/- CI
    ytop = Vvel_avg + CIfactor * Vvel_std;
    plot([x2 x2], [ybot ytop],       'Color',black,'LineWidth',2)
    plot([x2-dx x2+dx], [ybot ybot], 'Color',black,'LineWidth',2)
    plot([x2-dx x2+dx], [ytop ytop], 'Color',black,'LineWidth',2)
    hold off

    set(gca, 'Units', 'inches', 'YTick', YTickVals)
    set(gca, 'XTick', [1 2], 'XTickLabel', {'disp' 'vel'})
    set(gca, 'Box', 'on')
    ylim([0 maxPlotSpeed])
    xlim([0 3])
    ylabel('group SWS (m/s)')
    set(gca, 'Units', 'inches', 'Position', [0.5 0.5 1.4 2.7])

    kr = 2*pi*freqsToAnalyzeHz ./ phVel_avg .* rmin;        % kr values to find min freq in plot
    
    fmax = max(freqsToAnalyzeHz);                         % max freq based on size of error bars
    fbadidx = find(CIfactor*phVel_std>fracErrorBar*phVel_avg & kr>=krThreshold);
    if length(fbadidx) > 0
        fmax = freqsToAnalyzeHz(min(fbadidx) - 1);
    end

    fidx = GetPlotFreqIdx(kr, krThreshold, freqsToAnalyzeHz, fmax);

    capsize = 2;       % size of caps on black error bars
    patchwidth = 10;
    nn = length(fidx);
    patchx = NaN(1, 2*nn+2);
    patchy = NaN(1, 2*nn+2);

    for ii = 1:nn+1
        switch ii
            case 1
                patchx(ii)            = freqsToAnalyzeHz(fidx(ii)) - patchwidth/2;
                % if you are getting an error at this line, you may not
                % have a propogating shear wave, preventing further
                % analysis. This code assumes at least some phase 
                % velocities are in the valid range (and not omitted), if 
                % there aren't any valid points, it will error.
                % Check the planeDisp variable in FL#_phVel_gSWS_data.mat
                patchx(2*(nn+1)-ii+1) = freqsToAnalyzeHz(fidx(ii)) - patchwidth/2;
                patchy(ii)            = (1 + fracErrorBar) * phVel_avg(fidx(ii));
                patchy(2*(nn+1)-ii+1) = (1 - fracErrorBar) * phVel_avg(fidx(ii));
            case nn+1
                patchx(ii)            = freqsToAnalyzeHz(fidx(ii-1)) + patchwidth/2;
                patchx(nn+2)          = freqsToAnalyzeHz(fidx(ii-1)) + patchwidth/2;
                patchy(ii)            = (1 + fracErrorBar) * phVel_avg(fidx(ii-1));
                patchy(nn+2)          = (1 - fracErrorBar) * phVel_avg(fidx(ii-1));
            otherwise
                patchx(ii)            = freqsToAnalyzeHz(fidx(ii)) - patchwidth/2;
                patchx(2*(nn+1)-ii+1) = freqsToAnalyzeHz(fidx(ii)) - patchwidth/2;
                patchy(ii)            = (1 + fracErrorBar) / 2 * (phVel_avg(fidx(ii-1)) + phVel_avg(fidx(ii)));
                patchy(2*(nn+1)-ii+1) = (1 - fracErrorBar) / 2 * (phVel_avg(fidx(ii-1)) + phVel_avg(fidx(ii)));
        end
    end    

    hhh2 = subplot(1, 2, 2);
    hold on
    fill(patchx, patchy, [0.8 0.8 0.8], 'EdgeColor', 'none');

    for ii = 1:length(fidx)

        ybot = phVel_avg(fidx(ii)) - CIfactor * phVel_std(fidx(ii));        % +/- CI
        ytop = phVel_avg(fidx(ii)) + CIfactor * phVel_std(fidx(ii));
        plot([freqsToAnalyzeHz(fidx(ii))     freqsToAnalyzeHz(fidx(ii))    ], [ybot ytop], 'Color', black, 'LineWidth', 2)
        plot([freqsToAnalyzeHz(fidx(ii))-capsize freqsToAnalyzeHz(fidx(ii))+capsize], [ybot ybot], 'Color', black, 'LineWidth', 2)
        plot([freqsToAnalyzeHz(fidx(ii))-capsize freqsToAnalyzeHz(fidx(ii))+capsize], [ytop ytop], 'Color', black, 'LineWidth', 2)
    end
    hold off
    set(gca, 'Units', 'inches', 'YTick', YTickVals)
    set(gca, 'Box', 'on')
    ylim([0 maxPlotSpeed])
    xlim([0 maxplotfreq])
    ylabel('phase velocity (m/s)')
    xlabel('frequency (Hz)')
    title(phantomID)
    set(gca, 'Units', 'inches', 'Position', [2.4 0.5 5.4 2.7])

    savefig([saveDir '/gSWS_phVel_figure'])


        % make and save CSV file

    saveFile = [saveDir '/gSWS_phVel_data.txt'];

    fid = fopen(saveFile, 'w');
    fprintf(fid, '\n');
    fprintf(fid, 'phantom %s\n', phantomID);
    fprintf(fid, '\n');
    fprintf(fid, ',group SWS (m/s),  95%% CI,  30%% mean\n');
    fprintf(fid, 'displacement,%s,%s,%s\n', sprintf('%7.3f', Vdisp_avg),...
                                            sprintf('%7.3f', 1.95*Vdisp_std),...
                                            sprintf('%7.3f', 0.3*Vdisp_avg) );
    fprintf(fid, 'velocity,%s,%s,%s\n',     sprintf('%7.3f', Vvel_avg),...
                                            sprintf('%7.3f', 1.95*Vvel_std),...
                                            sprintf('%7.3f', 0.3*Vvel_avg) );
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    fprintf(fid, 'frequency (Hz),phase velocity (m/s),  95%% CI,  30%% mean\n');

    for ii = 1:length(fidx)
        fprintf(fid, '%s,%s,%s,%s\n', sprintf('%7.0f', freqsToAnalyzeHz(fidx(ii))),...
                                      sprintf('%7.3f', phVel_avg(fidx(ii))),...
                                      sprintf('%7.3f', 1.96*phVel_std(fidx(ii))),...
                                      sprintf('%7.3f', 0.30*phVel_avg(fidx(ii))) );
    end
    fclose(fid);

    errorFlag = 0;
    if (Vdisp_std*CIfactor > fracErrorBar*Vdisp_avg || Vvel_std*CIfactor > fracErrorBar*Vvel_avg ...
                || Vdisp_avg > par.maxValidGSWS || Vdisp_avg < par.minValidGSWS ...
                ||  Vvel_avg > par.maxValidGSWS ||  Vvel_avg < par.minValidGSWS)
        errorFlag=1;
    end
end

%==========================================================================

function [out_avg, out_std] = MeanStdOmitBadPtsOneCase_usingPositiveData(data, factor, maxSpeed)

    [npts, nfreqs] = size(data);
    out_avg = NaN(1, nfreqs);
    out_std = NaN(1, nfreqs);

    for ifreq = 1:nfreqs
        temp = data(:, ifreq);
        idx = find(temp>0 & temp<maxSpeed);
        temp_avg = mean(temp(idx));
        temp_std = std(temp(idx));

        out_avg(ifreq) = temp_avg;      % set outputs in case something goes wrong
        out_std(ifreq) = temp_std;      %   with goodidx below
        
        if temp_std > 0
            goodidx = find(abs(temp-temp_avg)<factor*temp_std & temp>0);        
        else
            goodidx = 1:npts;
        end
        out_avg(ifreq) = nanmean(temp(goodidx));
        out_std(ifreq) = nanstd(temp(goodidx));
    end
end

%==========================================================================

function [gSWS_avg, gSWS_std] = MeanStdOmitBadPtsOneCase_gSWS_usingPositiveData(data, factor, maxSpeed)
        
    idx = find(data>0 & data<maxSpeed);
    temp_avg = mean(data(idx));
    temp_std =  std(data(idx));

    gSWS_avg = temp_avg;      % set outputs in case something goes wrong
    gSWS_std = temp_std;      %   with goodidx below
        
    if temp_std > 0
        goodidx = find(abs(data-temp_avg)<factor*temp_std & data>0);        
    else
        goodidx = 1:length(data);
    end
    gSWS_avg = nanmean(data(goodidx));
    gSWS_std = nanstd(data(goodidx));
end

%==========================================================================

function freqIdx = GetPlotFreqIdx(kr, krThreshold, freqsToAnalyzeHz, fmax)
        
    freqIdx = find(kr>=krThreshold & freqsToAnalyzeHz<=fmax);
    if length(freqIdx) == 0
        fmax = max(freqsToAnalyzeHz);
        freqIdx = find(kr>=krThreshold & freqsToAnalyzeHz<=fmax);
    end
end

%==========================================================================
