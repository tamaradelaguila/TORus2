
for reject_on = [0 1]  %@ SET

user_settings;
 nfish =1;

clearvars -except nfish reject_on

VSDI = TORus('load',nfish);


% load waves to plot:
VSDroiTS= TORus('loadwave',nfish);
waves = VSDroiTS.filt306.data; %@ SET

Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
Lroi= [2 4 6 8 10 12 14];

setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET

pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_measurements/boxplot';

%% SELECT CASES
% cond_codes = unique(VSDI.condition(:,1));
% cond_codes=  cond_codes(~isnan(cond_codes));
% cond_codes= setdiff(cond_codes,0);

cond_codes =[300:302 ];

%% SELECT EXCLUDED

rejectidx = [];

if setting.manual_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
end

if setting.GSabsthres_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];
    
end

if setting.GSmethod_reject
    rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];
end

if setting.force_include
    rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
    
end

rejectidx = sort(unique(rejectidx));
    

%% SELECT CASES 
    sel_trials= [];
    
    for condi = cond_codes %to make sure that only the conditions of interest are computed)
        condtrials = makeCol(find(VSDI.condition(:,1)==condi));
        sel_trials  = [sel_trials; condtrials];
    end
    sel_trials= sort(sel_trials);
    
    if reject_on  %@ SET
        sel_trials = setdiff(sel_trials, rejectidx);
    end

    
    %% CALCULATE MEASURE FOR EACH TRIAL
    
        window.min = [-100 100];
        window.max = [0 600];
        window.movsum = 50;
        window.basel = [-100 0];
        window.slope = 50; %in ms
        
        method = 'movsum';
    

for nroi = 1:length(VSDI.roi.labels)
            for triali = makeRow(sel_trials)
                
                        wave = squeeze(waves(:, nroi,triali));
                        output = devo_peak2peak(wave, VSDI.timebase, window,[], method, 0);
                        
                        measures.peak2peak(nroi,triali) = output.p2p_value;
                        measures.peakminusbasel(nroi,triali) = output.peakminusbasel;
%                         frames.peaklat(rowi,coli,triali) = output.peaklat_ms;
%                         frames.p2plat(rowi,coli,triali) = output.p2plat_ms;
%                         frames.onset30_latency_ms(rowi,coli,triali) = output.onset30_latency_ms;
                        measures.onsetnoise_ms(nroi,triali) = output.onsetnoise_ms;
                        measures.noisethresh(nroi,triali) = output.noisethresh;
                        
                        
                                            
                        waveW = wave(output.peakidx(1):output.peakidx(2));
                        waveslope = diff(waveW);
                        meanslope = mean(waveslope);
                        
                        measures.meanslope(nroi,triali) = meanslope;

                        clear output wave waveW waveslope meanslope
                        
    
            end %triali
end %nroi
            
    
    
    
%% BOXPLOT FOR 'peakminusB'

       local_plotlim = [-.1 .4];
%        local_boxlim = ;

        %--------------------------------------
        % Right Hemisphere
        %--------------------------------------
     
        ploti = 1;%counter

    for nroi = makeRow(Rroi)
        measure = makeCol(squeeze(measures.peakminusbasel(nroi,sel_trials)));
        mA = VSDI.condition(sel_trials,4); %mA
        
        %
        subplot(3,3,ploti)
        boxplot(measure, mA, 'Colors', 'k')
        
%         ylim(local_boxlim);
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{nroi}])
        ploti = ploti+1;
        clear measure mA
        
    end %roi
    
    
%     subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
%     temp=  squeeze(mean(waves,2));
%     plot(VSDI.timebase,waves); hold on %plots the mean wave of all rois (GS from selected rois)
%     ylim(local_plotlim)
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    
    for nroi = Rroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
    axis image
    
    if reject_on
        sgtitle ([num2str(VSDI.ref), 'R: peakminusB ' '(cl)'])
    else
        sgtitle ([num2str(VSDI.ref), '_R: peakminusB'])
    end
    
    
    % Save
    if reject_on
        nameR = [num2str(VSDI.ref), 'boxplot peakminusB. Right(clean).jpg'];
    else
        nameR = [num2str(VSDI.ref), 'boxplot peakminusB. Right.jpg'];
    end
    saveas(gcf,fullfile(pathsave,nameR),'jpg')
    close
    
        %--------------------------------------
        % Left Hemisphere
        %--------------------------------------
     
        ploti = 1;%counter

    for nroi = makeRow(Lroi)
        measure = makeCol(squeeze(measures.peakminusbasel(nroi,sel_trials)));
        mA = VSDI.condition(sel_trials,4); %mA
        
        %
        subplot(3,3,ploti)
        boxplot(measure, mA, 'Colors', 'k')
        
%         ylim(local_boxlim);
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{nroi}])
        ploti = ploti+1;
        clear measure mA
        
    end %roi
    
    
%     subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
%     temp=  squeeze(mean(waves,2));
%     plot(VSDI.timebase,waves); hold on %plots the mean wave of all rois (GS from selected rois)
%     ylim(local_plotlim)
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    
    for nroi = Lroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
    axis image
    
    if reject_on
        sgtitle ([num2str(VSDI.ref), '_L: peakminusB ' '(cl)'])
    else
        sgtitle ([num2str(VSDI.ref), '_L: peakminusB'])
    end
    
    
    % Save
    if reject_on
        nameL = [num2str(VSDI.ref), 'boxplot peakminusB. Left(clean).jpg'];
    else
        nameL = [num2str(VSDI.ref), 'boxplot peakminusB. Left.jpg'];
    end
    saveas(gcf,fullfile(pathsave,nameL),'jpg')
    close
    
%% BOXPLOT FOR 'onsetnoise_ms'

       local_plotlim = [-.1 .4];
%        local_boxlim = ;

        %--------------------------------------
        % Right Hemisphere
        %--------------------------------------
     
        ploti = 1;%counter

    for nroi = makeRow(Rroi)
        measure = makeCol(squeeze(measures.onsetnoise_ms(nroi,sel_trials)));
        mA = VSDI.condition(sel_trials,4); %mA
        
        %
        subplot(3,3,ploti)
        boxplot(measure, mA, 'Colors', 'k')
        
%         ylim(local_boxlim);
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{nroi}])
        ploti = ploti+1;
        clear measure mA
        
    end %roi
    
    
%     subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
%     temp=  squeeze(mean(waves,2));
%     plot(VSDI.timebase,waves); hold on %plots the mean wave of all rois (GS from selected rois)
%     ylim(local_plotlim)
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    
    for nroi = Rroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
    axis image
    
    if reject_on
        sgtitle ([num2str(VSDI.ref), '_R: onsetnoise_m_s ' '(cl)'])
    else
        sgtitle ([num2str(VSDI.ref), '_R: onsetnoise_m_s'])
    end
    
    
    % Save
    if reject_on
        nameR = [num2str(VSDI.ref), 'boxplot onsetnoise_ms. Right(clean).jpg'];
    else
        nameR = [num2str(VSDI.ref), 'boxplot onsetnoise_ms. Right.jpg'];
    end
    saveas(gcf,fullfile(pathsave,nameR),'jpg')
    close
    
        %--------------------------------------
        % Left Hemisphere
        %--------------------------------------
     
        ploti = 1;%counter

    for nroi = makeRow(Lroi)
        measure = makeCol(squeeze(measures.onsetnoise_ms(nroi,sel_trials)));
        mA = VSDI.condition(sel_trials,4); %mA
        
        %
        subplot(3,3,ploti)
        boxplot(measure, mA, 'Colors', 'k')
        
%         ylim(local_boxlim);
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{nroi}])
        ploti = ploti+1;
        clear measure mA
        
    end %roi
    
    
%     subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
%     temp=  squeeze(mean(waves,2));
%     plot(VSDI.timebase,waves); hold on %plots the mean wave of all rois (GS from selected rois)
%     ylim(local_plotlim)
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    
    for nroi = Lroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
    axis image
    
    if reject_on
        sgtitle ([num2str(VSDI.ref), '_L: onsetnoise_m_s ' '(cl)'])
    else
        sgtitle ([num2str(VSDI.ref), '_L: onsetnoise_m_s'])
    end
    
    
    % Save
    if reject_on
        nameL = [num2str(VSDI.ref), 'boxplot onsetnoise_ms. Left(clean).jpg'];
    else
        nameL = [num2str(VSDI.ref), 'boxplot onsetnoise_ms. Left.jpg'];
    end
    saveas(gcf,fullfile(pathsave,nameL),'jpg')
    close
    

    %% BOXPLOT FOR 'meanslope'

       local_plotlim = [-.1 .4];
%        local_boxlim = ;

        %--------------------------------------
        % Right Hemisphere
        %--------------------------------------
     
        ploti = 1;%counter

    for nroi = makeRow(Rroi)
        measure = makeCol(squeeze(measures.meanslope(nroi,sel_trials)));
        mA = VSDI.condition(sel_trials,4); %mA
        
        %
        subplot(3,3,ploti)
        boxplot(measure, mA, 'Colors', 'k')
        
%         ylim(local_boxlim);
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{nroi}])
        ploti = ploti+1;
        clear measure mA
        
    end %roi
    
    
%     subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
%     temp=  squeeze(mean(waves,2));
%     plot(VSDI.timebase,waves); hold on %plots the mean wave of all rois (GS from selected rois)
%     ylim(local_plotlim)
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    
    for nroi = Rroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
    axis image
    
    if reject_on
        sgtitle ([num2str(VSDI.ref), '_R: meanslope ' '(cl)'])
    else
        sgtitle ([num2str(VSDI.ref), '_R: meanslope'])
    end
    
    
    % Save
    if reject_on
        nameR = [num2str(VSDI.ref), 'boxplot meanslope. Right(clean).jpg'];
    else
        nameR = [num2str(VSDI.ref), 'boxplot meanslope. Right.jpg'];
    end
    saveas(gcf,fullfile(pathsave,nameR),'jpg')
    close
    
        %--------------------------------------
        % Left Hemisphere
        %--------------------------------------
     
        ploti = 1;%counter

    for nroi = makeRow(Lroi)
        measure = makeCol(squeeze(measures.meanslope(nroi,sel_trials)));
        mA = VSDI.condition(sel_trials,4); %mA
        
        %
        subplot(3,3,ploti)
        boxplot(measure, mA, 'Colors', 'k')
        
%         ylim(local_boxlim);
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{nroi}])
        ploti = ploti+1;
        clear measure mA
        
    end %roi
    
    
%     subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
%     temp=  squeeze(mean(waves,2));
%     plot(VSDI.timebase,waves); hold on %plots the mean wave of all rois (GS from selected rois)
%     ylim(local_plotlim)
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    
    for nroi = Lroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
    axis image
    
    if reject_on
        sgtitle ([num2str(VSDI.ref), '_L: meanslope ' '(cl)'])
    else
        sgtitle ([num2str(VSDI.ref), '_L: meanslope'])
    end
    
    
    % Save
    if reject_on
        nameL = [num2str(VSDI.ref), 'boxplot meanslope. Left(clean).jpg'];
    else
        nameL = [num2str(VSDI.ref), 'boxplot meanslope. Left.jpg'];
    end
    saveas(gcf,fullfile(pathsave,nameL),'jpg')
    close
    
    
    
end 