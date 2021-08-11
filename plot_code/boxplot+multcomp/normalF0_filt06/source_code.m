%% TEST BOXPLOT + STATS
clear
close all

nfish = 4;
reject_on = 1;
user_settings;

%load structure + wave
VSDI = TORus('load',nfish);
VSDroiTS= TORus('loadwave',nfish);

% roikind = 'circle';
roikind = 'anat';


cond_codes =[100:102];

% PARAMETERS

window.min = [-100 100];
window.max = [0 600];
window.movsum = 50;
window.basel = [-100 0];

method = 'movsum';

%%
Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
Lroi= [2 4 6 8 10 12 14];

switch roikind
    case 'circle'
        waves = VSDroiTS.circ_filt306.data; %@ SET
    case 'anat'
        waves = VSDroiTS.filt306.data; %@ SET
end

setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET



% SELECT EXCLUDED

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


% SELECT CASES
sel_trials= [];

for condi = cond_codes %to make sure that only the conditions of interest are computed)
    condtrials = makeCol(find(VSDI.condition(:,1)==condi));
    sel_trials  = [sel_trials; condtrials];
end
sel_trials= sort(sel_trials);

if reject_on  %@ SET
    sel_trials = setdiff(sel_trials, rejectidx);
end


% CALCULATE MEASURE FOR EACH TRIAL

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
% 
% subplot(1,5,1)
% idx = find(VSDI.condition(:,1) == 400); idx = intersect(idx,sel_trials );
%        plot(VSDI.timebase, squeeze(mean(waves(:, :,idx),3))); ylim([-0.1 0.6]); xline(0)
%       
%        subplot(1,5,2)
% idx = find(VSDI.condition(:,1) == 401); idx = intersect(idx,sel_trials );
%        plot(VSDI.timebase, squeeze(mean(waves(:, :,idx),3))); ylim([-0.1 0.6]); xline(0)
% 
%               subplot(1,5,3)
% idx = find(VSDI.condition(:,1) == 402); idx = intersect(idx,sel_trials );
%        plot(VSDI.timebase, squeeze(mean(waves(:, :,idx),3))); ylim([-0.1 0.6]); xline(0)
% 
%                      subplot(1,5,4)
% idx = find(VSDI.condition(:,1) == 403); idx = intersect(idx,sel_trials );
%        plot(VSDI.timebase, squeeze(mean(waves(:, :,idx),3))); ylim([-0.1 0.6]); xline(0)
% 
%                      subplot(1,5,5)
% idx = find(VSDI.condition(:,1) == 404); idx = intersect(idx,sel_trials );
%        plot(VSDI.timebase, squeeze(mean(waves(:, :,idx),3))); ylim([-0.1 0.6]); xline(0)

%% BOXPLOT + STAT
close all
roiname = 'dm2_R';
% result = 'peakminusbasel';
% result ='onsetnoise_ms';
result = 'meanslope';

boxplot_ylim = [-.01 .04];

% close all


%%
nroi = name2idx(roiname, VSDI.roi.labels);

measure = [];
switch result
    case 'peakminusbasel'
        measure = makeCol(squeeze(measures.peakminusbasel(nroi,sel_trials)));
    case 'onsetnoise_ms'
        measure = makeCol(squeeze(measures.onsetnoise_ms(nroi,sel_trials)));
    case 'meanslope'
        measure = makeCol(squeeze(measures.meanslope(nroi,sel_trials)));
end
mA= [];
mA = VSDI.condition(sel_trials,4); %mA

% BOXPLOT

subplot(1,3,1)
boxplot(measure, mA, 'Colors', 'k')
% ylim(boxplot_ylim);

xlabel('mA'); ylabel('% /Delta F');

title([ VSDI.roi.labels{nroi}, ':', result])


% subplot(1,3,2) %plot the R-hemisph in the last empty plot (to have a visual guide)
% temp=  squeeze(squeeze(mean(waves(:,nroi,sel_trials),2)));
% plot(VSDI.timebase,temp'); hold on %plots the mean wave of all rois (GS from selected rois)
% %     ylim(local_plotlim)


switch roikind
    case 'circle'
        axH= subplot(1,2,2); %plot the R-hemisph in the last empty plot (to have a visual guide
        
        roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(Rroi,:), VSDI.roi.circle.R, axH);
        title(['#' num2str(VSDI.ref)])
        axis image
        
    case 'anat'
        subplot(1,2,2);
        imagesc(VSDI.crop.preview); colormap('bone'); hold on
        roicolors= roi_colors();
        
            coord = VSDI.roi.manual_poly{nroi};
            fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
            plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
            title(['#' num2str(VSDI.ref)])

            axis image
        
end



%STATS

[p, table, stats] = anova1(measure,mA, 'on')
figure
[c, m , h , gnames] = multcompare(stats)
title([ num2str(VSDI.ref), VSDI.roi.labels{nroi}])
% end stats

% close all

% %% MAKE long format
% 
% labels = unique(mA)
% for ii  =  makeRow(labels) %to make sure that only the conditions of interest are computed)
% tidx = find(mA == ii);
% temp =  ['mA' num2str(ii)]; temp = erase(temp, '.')
%     spss_peak.(temp)= measure(tidx);
% end
% 
% 

% %% FOR FISH 6 (#210412): instead of mA, group by cond_code
% nroi = name2idx(roiname, VSDI.roi.labels);
% 
% measure = [];
% switch result
%     case 'peakminusbasel'
%         measure = makeCol(squeeze(measures.peakminusbasel(nroi,sel_trials)));
%     case 'onsetnoise_ms'
%         measure = makeCol(squeeze(measures.onsetnoise_ms(nroi,sel_trials)));
%     case 'meanslope'
%         measure = makeCol(squeeze(measures.meanslope(nroi,sel_trials)));
% end
% code= [];
% code = VSDI.condition(sel_trials,1); %mA
% 
% % BOXPLOT
% 
% subplot(1,3,1)
% boxplot(measure, code, 'Colors', 'k')
% % ylim(boxplot_ylim);
% 
% xlabel('cond'); ylabel('% /Delta F');
% 
% title([ VSDI.roi.labels{nroi}, ':', result])
% 
% 
% % subplot(1,3,2) %plot the R-hemisph in the last empty plot (to have a visual guide)
% % temp=  squeeze(squeeze(mean(waves(:,nroi,sel_trials),2)));
% % plot(VSDI.timebase,temp'); hold on %plots the mean wave of all rois (GS from selected rois)
% % %     ylim(local_plotlim)
% 
% 
% switch roikind
%     case 'circle'
%         axH= subplot(1,2,2); %plot the R-hemisph in the last empty plot (to have a visual guide
%         
%         roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(Rroi,:), VSDI.roi.circle.R, axH);
%         title(['#' num2str(VSDI.ref)])
%         axis image
%         
%     case 'anat'
%         subplot(1,2,2);
%         imagesc(VSDI.crop.preview); colormap('bone'); hold on
%         roicolors= roi_colors();
%         
%             coord = VSDI.roi.manual_poly{nroi};
%             fill(coord(:,1), coord(:,2), roicolors(nroi,:),'FaceAlpha',0.5); hold on
%             plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
%             title(['#' num2str(VSDI.ref)])
% 
%             axis image
%         
% end
% 
% 
% 
% %STATS
% 
% [p, table, stats] = anova1(measure,code, 'on')
% figure
% [c, m , h , gnames] = multcompare(stats)
% title([ num2str(VSDI.ref), VSDI.roi.labels{nroi}])
% % end stats
% 
% % close all

