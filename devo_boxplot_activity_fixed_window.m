%% 1. COMPARE DIFFERENT PARAMETERS 

clear
user_settings;
nfish =4;

VSDI = TORus('load',nfish);
% spike = TORus('loadspike', nfish);

% load waves to plot: 
VSDroiTS =TORus('loadwave',nfish);
timeseries = VSDroiTS.filt306.data;  %@ SET
Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
Lroi= [2 4 6 8 10 12 14];

manual_reject =1; %@ SET

%% SELECT CASES: OPTION 1 (manually set specific conditions)
c1 = 100; c2 = 103;% @SET find extreme codes to select from VSDI.list.code
sel_trials = find(VSDI.condition(:,1)>=c1 & VSDI.condition(:,1)<=c2);
    if manual_reject
       sel_trials= setdiff(sel_trials, VSDI.reject.manual);
    end
temp = num2str(c1); cond_def= [temp(1),'-',temp(2)]; clear temp %to get the code of the stimulus

%% 1.1 BUILD ACTIVITY MATRIX (trials - roi)
% Measure: mean activity in timewindow
% wind_ms = [202 502]; %@ SET
windows = [0 100; 0 200; 100 200; 100 300] ; %@SET

for wi = 1:size(windows,1) %will loop through rows
    
%1. WINDOW
% wind_ms = [202 502]; %@ SET
wind_ms = windows(wi,:);
 
wind_idx = find_closest_timeidx(wind_ms, VSDI.timebase);
wind_actualms =  VSDI.timebase(wind_idx);

roi_mean_act= mean (timeseries(wind_idx(1):wind_idx(2), :,:)); %roi_mean_act: roi x trials
roi_mean_act = squeeze(roi_mean_act); %to avoid having a first unique dimension
%roi_mean_act dimensions are : roixtrial

% Measure: calculate slope in that time window
% 
% temp = timeseries(wind_idx(1):wind_idx(2), : , :); 
% temp2 = diff(temp,1); 
% roi_mean_slope = squeeze(mean(temp2,1)); %dimensions roi_mean_slope: roixtrials
% 
% for triali = VSDI.nonanidx
%     for roi = 1:14
%         t = VSDI.timebase(wind_idx(1):wind_idx(2));
%         y = timeseries(wind_idx(1):wind_idx(2),roi, triali);
% 
%         p = polyfit(t,y,1);
%         y_est =polyval(p,t);
%         plot(t,y); hold on
%         plot(t,y_est)
%         
%         slope(roi,triali)=  atan(p(1)); %da el ángulo en radianes
%         slope(roi,triali)=  atand(p(1)); %da el ángulo en radianes
% 
%     end
% end
% clear temp1 temp2


%% 1.2.1 BOXPLOT FROM MEAN_ ACTIVITY

savepath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/'; %@ SET
act_lim = [-0.1 0.4];

% RIGHT HEMISPH:
figure
ploti = 1;
for roi = Rroi

activ = roi_mean_act(roi,sel_trials);
mA = VSDI.condition(sel_trials,4);

% condition copy: to fix that the tone has no mA (we need a number for the
% boxplot
temp_cond = VSDI.condition(sel_trials,2);
mA(isnan(temp_cond)) = 600; %to indicate the tone - check if it can be done 

% 
subplot(3,3,ploti)
boxplot(activ, mA, 'Colors', 'k')
ylim(act_lim)
xlabel('mA'); ylabel('% /Delta F'); 
title([ VSDI.roi.labels{roi} ':' cond_def])
ploti = ploti+1;
clear act mA 
end

subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
temp=  squeeze(mean(timeseries,2));
plot(VSDI.timebase,temp); hold on %plots the mean wave of all rois (GS from selected rois)
ylim(act_lim)
down= act_lim(1);up= act_lim(2); 
patch ([wind_ms(1) wind_ms(2), wind_ms(2), wind_ms(1)],[ down down up up],'k','FaceAlpha',.3, 'LineStyle', 'none')


subplot(3,3,9) %plot the R-hemisph in the last empty plot (to have a visual guide)
    imagesc(VSDI.crop.preview); colormap('bone'); hold on
    roicolors= roi_colors();
    for nroi = Rroi
        coord = VSDI.roi.manual_poly{nroi,1};
        fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end
if manual_reject
        sgtitle ([num2str(VSDI.ref), 'R_H: mean in window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms. (cl)'])
    else 
        sgtitle ([num2str(VSDI.ref), 'R_H: mean in window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms'])
    end


% Save 
    if manual_reject
        nameR = [num2str(VSDI.ref), 'cond', cond_def,'mean_window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - Right(cl)'];
    else 
        nameR = [num2str(VSDI.ref), 'cond', cond_def,'mean_window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - Right'];
    end
    saveas(gcf,fullfile(savepath,nameR),'jpg')
    close
    
% LEFT HEMISPH:
figure
ploti = 1;
for roi = Lroi

activ = roi_mean_act(roi,sel_trials);
mA = VSDI.condition(sel_trials,4);
% 
% we fix that the tone has no mA (we need a number for the
% boxplot
temp_cond = VSDI.condition(sel_trials,2);
mA(isnan(temp_cond)) = 600; 
subplot(3,3,ploti)
boxplot(activ, mA, 'Colors', 'k')
ylim(act_lim)
xlabel('mA'); ylabel('% /Delta F'); 
title([ VSDI.roi.labels{roi} ':' cond_def])
ploti = ploti+1;
    if manual_reject
        sgtitle ([num2str(VSDI.ref), 'L_H: mean in window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms. (cl)'])
    else 
        sgtitle ([num2str(VSDI.ref), 'L_H: mean in window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms'])
    end
clear act mA 
end

subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
temp=  squeeze(mean(timeseries,2));
plot(VSDI.timebase,temp); hold on %plots the mean wave of all rois (GS from selected rois)
ylim(act_lim)
down= act_lim(1);up= act_lim(2); 
patch ([wind_ms(1) wind_ms(2), wind_ms(2), wind_ms(1)],[ down down up up],'k','FaceAlpha',.3, 'LineStyle', 'none')

subplot(3,3,8) %plot the R-hemisph in the last empty plot (to have a visual guide)
temp=  squeeze(mean(timeseries,2));
plot(VSDI.timebase,temp); hold on %plots the mean wave of all rois (GS from selected rois)
ylim(act_lim)
down= act_lim(1);up= act_lim(2); 
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
    if manual_reject
        nameL = [num2str(VSDI.ref), 'cond', cond_def,'mean_window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - Left(cl)'];
    else 
        nameL = [num2str(VSDI.ref), 'cond', cond_def,'mean_window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms - Left'];
    end

    saveas(gcf,fullfile(savepath,nameL),'jpg')
    close
    
    clear nameR nameL savepathR savepathL 
    
end % wi
% 
% % to perform multiple comparisons: choose one roi
% %% 1.3 MULTIPLE COMPARISONS  statistics for each roi (check the mult comparisons table and then press any key in the command window to advance)
% 
% for roi = [3 4 7 8] %roi2plotidx
% activ = roi_mean_act(roi,sel_trials);
% mA = VSDI.condition(sel_trials,4);
% 
% [p, table, stats] = anova1(activ,mA, 'on')
% [c, m , h , gnames] = multcompare(stats)
% title(VSDI.roi.labels{roi})
% pause
% close all
% 
% end
% 
% %% 2. COMPARE ROI for each Intensity
% 
% %% 2.1 BUILD ACTIVITY MATRIX (trials - roi)
% % Measure: mean activity in timewindow
% wind_ms = [202 502]; %@ SET
% wind_idx = find_closest_timeidx(wind_ms, VSDI.timebase);
% wind_actualms =  VSDI.timebase(wind_idx);
% 
% roi_mean_act= mean (timeseries(wind_idx(1):wind_idx(2), :,:)); % roi x trials
% roi_mean_act = squeeze(roi_mean_act); %to avoid having a first unique dimension
% %roi_mean_act dimensions are : roixtrial
% 
% 
% % BUILD LONG FORMAT 
% ntrial = length(VSDI.list);
% i = 1;
% for nroi = 1:length(VSDI.roi.labels)
%     long(i:i+ntrial-1,1) = nroi;% roi index
%     long(i:i+ntrial-1,2) = roi_mean_act(nroi,:); % mean activity of the roi for each trial
%     long(i:i+ntrial-1,3) = VSDI.condition(:,1); %code (to select trials)
%     long(i:i+ntrial-1,4) = VSDI.condition(:,4); %intensity
%     i = i+ntrial;
%     
% end
% 
% %% 2.2 BOXPLOT 
% % Arrays for boxplot
% figure
% ploti = 1;
% for cond = 100:103
% 
% sel_idx = find(long(:,3)==cond);
% cond_def = [num2str(cond) ':' num2str(long(sel_idx(1),4)) 'mA']; %to get the code of the stimulus
% 
% activ = long(sel_idx,2);
% roi = long(sel_idx, 1);
% 
% %turn into string labels to input in the 
% for ii = 1:length(roi)
%     roi_idx = roi(ii);
%     roilab{ii} = VSDI.roi.labels{roi_idx};
% end
% 
% % 
% subplot(4,1,ploti)
% boxplot(activ, roilab, 'Colors', 'k')
% ylim([-0.1 0.7])
% ylabel('mean % /Delta F');
% title(cond_def)
% xtickangle(45)
% xticklabels(VSDI.roi.labels);
% % title([ VSDI.roi.labels{roi} ':' cond_def])
% ploti = ploti+1;
% sgtitle ([num2str(VSDI.ref), 'mean in window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms'])
% clear roilab roi_idx roi activ
% end
% 
% % Save 
% %     name = [num2str(VSDI.ref), 'cond', cond_def,'mean_window:', num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'ms'];
% %     savepath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/mean_act_wind200to500ms';
% %     saveas(gcf,fullfile(savepath,name),'jpg')
% %     close 
% 
% %% 2.3 %% 1.3 MULTIPLE COMPARISONS  statistics for each condition (check the mult comparisons table and then press any key in the command window to advance) 
% for cond = 100:103
% 
% sel_idx = find(long(:,3)==cond);
% cond_def = [num2str(cond) ':' num2str(long(sel_idx(1),4)) 'mA']; %to get the code of the stimulus
% 
% activ = long(sel_idx,2);
% roi = long(sel_idx, 1);
% 
% %turn into string labels to input in the 
% for ii = 1:length(roi)
%     roi_idx = roi(ii);
%     roilab{ii} = VSDI.roi.labels{roi_idx};
% end
% 
% 
% [p, table, stats] = anova1(activ,roilab, 'on')
% [c, m , h , gnames] = multcompare(stats)
% title([num2str(VSDI.ref) '-' num2str(cond_def)])
% pause
% close all
%  
% clear activ roi roi_idx roilab
% end
