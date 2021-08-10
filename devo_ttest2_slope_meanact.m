clear

user_settings;
for nfish = 6

clearvars -except nfish

VSDI = TORus('load',nfish);
% spike = TORus('loadspike', nfish);

% load waves to plot:
VSDroiTS =TORus('loadwave',nfish);
timeseries = VSDroiTS.filt306.data;  %@ SET
Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
Lroi= [2 4 6 8 10 12 14];

setting.reject_on = 1;  %@ SET
setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET

pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/slope_mean_pmap';

%% SELECT CASES
cond_codes = unique(VSDI.condition(:,1));
cond_codes=  cond_codes(~isnan(cond_codes));
cond_codes= setdiff(cond_codes,0);

% cond_codes =[2000 2001 2002 2003];

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

% %% TTEST ROI-TO-ROI COMPARISONS 
% % Measure: mean activity in timewindow
% % wind_ms = [202 502]; %@ SET
% windows = [60 160] ; %@SET
% 
% for wi= 1:size(windows,1)
%     wind_ms = windows(wi,:);
%     % wind_ms = windows(wi,:);
%     
%     wind_idx = find_closest_timeidx(wind_ms, VSDI.timebase);
%     wind_actualms =  VSDI.timebase(wind_idx);
%     
%     idx_range = wind_idx(1):wind_idx(end);
%     
%     waves = timeseries(idx_range,:,:);
%     wavesSlope = diff(waves, [], 1);
%     
%     roi_mean_slope= squeeze(mean(wavesSlope));%for each trial;
%     
%     
%     %% BUILD MATRIX FOR MEAN ACTIVITY
%     
%     roi_mean_act= mean (timeseries(wind_idx(1):wind_idx(2), :,:)); %roi_mean_act: roi x trials
%     roi_mean_act = squeeze(roi_mean_act); %to avoid having a first unique dimension
%     
%     %% FOR EACH CONDITION, CALCULATE MATRIX AND PLOT
%     for ci =1:numel(cond_codes) %for each condition
%         
%         sel_trials = find(VSDI.condition(:,1) == cond_codes(ci));
%         if setting.reject_on
%             sel_trials= setdiff(sel_trials, rejectidx);
%         end
%         
%         cond_def = num2str(cond_codes(ci));
%         
%         % CALCULATE MATRIX OF P-VALUES: roi-roi matrix of relevant differences in the meassure
%         nroi = length(VSDI.roi.labels);
%         
%         for roi1 = 1:nroi
%             for roi2 = roi1:nroi %only for one half
%                 
%                 % SLOPE
%                 s1 = roi_mean_slope(roi1,sel_trials); %measure for roi 1
%                 s2 = roi_mean_slope(roi2,sel_trials); %measure for roi 1
%                 
%                 h_slope(roi2,roi1) = ttest2(s1,s2, 'alpha', 0.05);
%                 
%                 % MEAN
%                 m1 = roi_mean_act(roi1,sel_trials); %measure for roi 1
%                 m2 = roi_mean_act(roi2,sel_trials); %measure for roi 1
%                 
%                 h_mean(roi2,roi1) = ttest2(m1,m2, 'alpha', 0.05);
%                 
%             end
%         end
%         
%         % PLOT SLOPE
%         figure
%         for roiplot = 1:14
%         subplot (4,4,ploti)
%         imagesc(h_slope(roiplot,:,:));
%         % colormap = polarmap();
%         
%                 
%         xticks(1:nroi)
%         xticklabels (roilabels); xtickangle(90)
%         yticks(1:nroi)
%         yticklabels (roilabels)
%         set(gca, 'clim', [0 1]);
%         
%         end
%         sgtitle ('mean slope')
% 
%         
%         % PLOT MEAN
%         subplot (1,2,2)
%         imagesc(h_mean);
%         % colormap = polarmap();
%         title ('mean activity')
%         % set(gca,'clim', climit);
%         set(gca, 'clim', [0 1])
%         
%         xticks(1:nroi)
%         xticklabels (roilabels); xtickangle(90)
%         yticks(1:nroi)
%         yticklabels (roilabels)
%         
%         sgtitle([num2str(VSDI.ref), '.w=' ,num2str(wind_actualms(1)), 'to', num2str(wind_actualms(2)), ' (cond' cond_def, '). dif p<0.05']);
%         
%         name = [num2str(VSDI.ref),'(w=',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)), 'ms: cond=',num2str(cond_codes(ci))];
%         saveas(gcf, fullfile(pathsave, name), 'jpg')
%         
%     end %loop ci
% end %loop wi

%% TTEST ROI-CONDITION COMPARISON
%% TTEST ROI-CONDITION COMPARISON

windows = [60 160] ; %@SET
pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/ttest2_conditions';

for wi= 1:size(windows,1)
    wind_ms = windows(wi,:);
    % wind_ms = windows(wi,:);
    
    wind_idx = find_closest_timeidx(wind_ms, VSDI.timebase);
    wind_actualms =  VSDI.timebase(wind_idx);
    
    idx_range = wind_idx(1):wind_idx(end);
    
    waves = timeseries(idx_range,:,:);
    wavesSlope = diff(waves, [], 1);
    
    roi_mean_slope= squeeze(mean(wavesSlope));%for each trial;
    
    
    %% BUILD MATRIX FOR MEAN ACTIVITY
    
    roi_mean_act= mean (timeseries(wind_idx(1):wind_idx(2), :,:)); %roi_mean_act: roi x trials
    roi_mean_act = squeeze(roi_mean_act); %to avoid having a first unique dimension
    
    %% FOR EACH CONDITION, CALCULATE MATRIX AND PLOT
        
        
        
        % CALCULATE MATRIX OF P-VALUES: roi-roi matrix of relevant differences in the meassure
        nroi = length(VSDI.roi.labels);
        
        for roi = 1:nroi
            
            for c1 = 1:numel(cond_codes)
                
                sel_trials1 = find(VSDI.condition(:,1) == cond_codes(c1));
                if setting.reject_on
                    sel_trials1= setdiff(sel_trials1, rejectidx);
                end
                
                for c2 = c1:numel(cond_codes)
                    
                    sel_trials2 = find(VSDI.condition(:,1) == cond_codes(c2));
                    if setting.reject_on
                        sel_trials2= setdiff(sel_trials2, rejectidx);
                    end
                    
                    
                    % SLOPE
                    s1 = roi_mean_slope(roi,sel_trials1); %measure for roi 1
                    s2 = roi_mean_slope(roi,sel_trials2); %measure for roi 1
                    
                    h_slope(roi, c1,c2) = ttest2(s1,s2, 'alpha', 0.05);
                    
                    % MEAN
                    m1 = roi_mean_act(roi,sel_trials1); %measure for roi 1
                    m2 = roi_mean_act(roi,sel_trials2); %measure for roi 1
                    
                    h_mean(roi,c1, c2) = ttest2(m1,m2, 'alpha', 0.05);
                    
                end
            end
        end
        
        % BUILD CONDITIONS LABELS
        
        ncond = length(cond_codes);
        for ii=1:ncond
           cond_label{ii}=num2str(cond_codes(ii)); 
        end

        
        % PLOT SLOPE
        figure
        if setting.reject_on
            
            sgtitle([num2str(VSDI.ref), 'ttest2 SLOPE.w=' ,num2str(wind_actualms(1)), 'to', num2str(wind_actualms(2)), 'ms.dif p<0.05 (clean)']);
            name1 = [num2str(VSDI.ref),'ttest2conditions_slope.w=',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)), 'ms_clean.jpg'];
        else
            sgtitle([num2str(VSDI.ref), 'ttest2 SLOPE.w=' ,num2str(wind_actualms(1)), 'to', num2str(wind_actualms(2)), 'ms.dif p<0.05']);
            name1 = [num2str(VSDI.ref),'ttest2conditions_slope.w=',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)), 'ms.jpg'];
        end
                
        for roiplot = 1:14
        subplot (4,4,roiplot)
        imagesc(squeeze(h_slope(roiplot,:,:)));
        % colormap = polarmap();
        % set(gca,'clim', climit);
        title(VSDI.roi.labels{roiplot})
        
        xticks(1:ncond)
        xticklabels (cond_codes); xtickangle(90)
        yticks(1:ncond)
        yticklabels (cond_codes)
        set(gca,'fontsize',6)

        set(gca, 'clim', [0 1]);
        end
         set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 6]); %

        saveas(gcf, fullfile(pathsave, name1), 'jpg')
        close

        
        % PLOT MEAN
        figure
        if setting.reject_on
            sgtitle([num2str(VSDI.ref), 'ttest2 ACT.w=' ,num2str(wind_actualms(1)), 'to', num2str(wind_actualms(2)), 'ms.dif p<0.05 (clean)']);
            name2 = [num2str(VSDI.ref),'ttest2conditions_Act.w=',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)), 'ms_clean.jpg'];

        else
            sgtitle([num2str(VSDI.ref), 'ttest2 ACT.w=' ,num2str(wind_actualms(1)), 'to', num2str(wind_actualms(2)), 'ms.dif p<0.05']);
            name2 = [num2str(VSDI.ref),'ttest2conditions_Act.w=',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)), 'ms.jpg'];

        end
        
        
        for roiplot = 1:14
            subplot (4,4,roiplot)
            imagesc(squeeze(h_mean(roiplot,:,:)));
            % colormap = polarmap();
            % set(gca,'clim', climit);
            title(VSDI.roi.labels{roiplot})
            
            xticks(1:ncond)
            xticklabels (cond_codes); xtickangle(90)
            yticks(1:ncond)
            yticklabels (cond_codes)
            set(gca,'fontsize',6)

            set(gca, 'clim', [0 1]);
        end
        
         set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 6 6]); %

        saveas(gcf, fullfile(pathsave, name2), 'jpg')
        close
end %loop wi
end % nfish