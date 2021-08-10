%% BOXPLOT SLOPE IN FIXED WINDOW
% Conditions to be include in each has to be set for each condition
% manually, if we want them grouped by conditions' blocks


clear
user_settings;
nfish = 6; % 2 3 4 8

VSDI = TORus('load',nfish);
% spike = TORus('loadspike', nfish);

% load waves to plot:
VSDroiTS =TORus('loadwave',nfish);
timeseries = VSDroiTS.filt306.data;  %@ SET
Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
Lroi= [2 4 6 8 10 12 14];

setting.reject_on = 0;  %@ SET
setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 0; %@ SET

%% SELECT CASES: OPTION 1 (manually set specific conditions)
cond_codes = sort(unique(VSDI.condition(:,1))); 
cond_codes = cond_codes(~isnan(cond_codes))

c1 = 4000; c2 = 4002;% @SET find extreme codes to select from VSDI.list.code
temp = num2str(c1); cond_def= [temp(1),'-',temp(2)]; clear temp %to get the code of the stimulus

sel_trials = find(VSDI.condition(:,1)>=c1 & VSDI.condition(:,1)<=c2);


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

if setting.reject_on  %@ SET
    
sel_trials = setdiff(sel_trials, rejectidx);
    
end

%% 1.1 BUILD ACTIVITY MATRIX (trials - roi)
% Measure: mean activity in timewindow
% wind_ms = [202 502]; %@ SET
windows = [60 160; 60 210] ; %@SET

for wi = 1:size(windows,1) %will loop through rows
    
    %1. WINDOW
    % wind_ms = [202 502]; %@ SET
    wind_ms = windows(wi,:);
    
    wind_idx = find_closest_timeidx(wind_ms, VSDI.timebase);
    wind_actualms =  VSDI.timebase(wind_idx);
    
    idx_range = wind_idx(1):wind_idx(end);
    
    waves = timeseries(idx_range,:,:);
    wavesSlope = diff(waves, [], 1);
    
    roi_mean_slope= squeeze(mean(wavesSlope));%for each trial;
    
    
    %% 1.2.1 BOXPLOT FROM MEAN_SLOPE
    
    savepath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/boxplot_03'; %@ SET
    
    %3. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
    % temp_include = roi_mean_slope([Rroi, Lroi],sel_trials);
    % maxval= max(abs(temp_include(:))); %get extreme values from included waves and trials
    % slope_lim = [-maxval maxval];
    slope_lim = [-.005 0.015];
    plot_lim = [-.1 .4];
    BVmap = colormap_loadBV();
    
    % RIGHT HEMISPH:
    figure
    ploti = 1;%counter
    for roi = Rroi
        
        activ = roi_mean_slope(roi,sel_trials);
        mA = VSDI.condition(sel_trials,4);
        
        % IF NEEDED
        % condition copy: to fix that the tone has no mA (we need a number for the
        % boxplot
        temp_cond = VSDI.condition(sel_trials,2);
        mA(isnan(temp_cond)) = 600; %to indicate the tone - check if it can be done
        
        %
        subplot(3,3,ploti)
        boxplot(activ, mA, 'Colors', 'k')
        ylim(slope_lim)
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{roi} ':' cond_def])
        ploti = ploti+1;
        clear act mA
    end
    
    subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
    temp=  squeeze(mean(timeseries,2));
    plot(VSDI.timebase,temp); hold on %plots the mean wave of all rois (GS from selected rois)
    ylim(plot_lim)
    down= plot_lim(1);up= plot_lim(2);
    patch ([wind_ms(1) wind_ms(2), wind_ms(2), wind_ms(1)],[ down down up up],'k','FaceAlpha',.3, 'LineStyle', 'none')
    
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    for nroi = Rroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
    if setting.reject_on
        sgtitle ([num2str(VSDI.ref), 'R_H: mean slope in w=', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms. (cl)'])
    else
        sgtitle ([num2str(VSDI.ref), 'R_H: mean slope in w=', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms'])
    end
    
    
    % Save
    if setting.reject_on
        nameR = [num2str(VSDI.ref), 'boxplot slope. cond', cond_def,':', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - cond', num2str(c1) ,'-Right(clean).jpg'];
    else
        nameR = [num2str(VSDI.ref), 'boxplot slope. cond', cond_def,':', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - cond', num2str(c1) , 'Right.jpg'];
    end
    saveas(gcf,fullfile(savepath,nameR),'jpg')
    close
    
    % LEFT HEMISPH:
    figure
    ploti = 1;
    for roi = Lroi
        
        activ = roi_mean_slope(roi,sel_trials);
        mA = VSDI.condition(sel_trials,4);
        %
        % we fix that the tone has no mA (we need a number for the
        % boxplot
        temp_cond = VSDI.condition(sel_trials,2);
        mA(isnan(temp_cond)) = 600;
        subplot(3,3,ploti)
        boxplot(activ, mA, 'Colors', 'k')
        ylim(slope_lim)
        xlabel('mA'); ylabel('% /Delta F');
        title([ VSDI.roi.labels{roi} ':' cond_def])
        ploti = ploti+1;
        if setting.reject_on
            sgtitle ([num2str(VSDI.ref), 'L_H: mean slope in w=', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms. (cl)'])
        else
            sgtitle ([num2str(VSDI.ref), 'L_H: mean slope in w=', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms'])
        end
        clear activ mA
    end
    
    
    subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
    temp=  squeeze(mean(timeseries,2));
    plot(VSDI.timebase,temp); hold on %plots the mean wave of all rois (GS from selected rois)
    ylim(plot_lim)
    down= plot_lim(1);up= plot_lim(2);
    patch ([wind_ms(1) wind_ms(2), wind_ms(2), wind_ms(1)],[ down down up up],'k','FaceAlpha',.3, 'LineStyle', 'none')
    
    
    subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    
    for nroi = Lroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
        
    end
    
    
    % Save L-hemisph
    if setting.reject_on  %@ SET
        
        nameL = [num2str(VSDI.ref), 'boxplot slope. cond', cond_def,':', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - cond', num2str(c1) , 'Left(clean).jpg'];
    else
        nameL = [num2str(VSDI.ref), 'boxplot slope. cond', cond_def,':', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - cond', num2str(c1) , 'Left.jpg'];
    end
    
    saveas(gcf,fullfile(savepath,nameL),'jpg')
    close
    
    clear nameR nameL savepathR savepathL
    
end % wi
