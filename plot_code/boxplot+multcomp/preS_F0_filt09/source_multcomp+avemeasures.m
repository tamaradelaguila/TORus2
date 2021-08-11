%% TEST BOXPLOT + STATS

% SETTINGS
clear
% close all

nfish = 11;
reject_on = 1;
user_settings;

%load structure + wave
VSDI = TORus('load',nfish);
VSDroiTS= TORus('loadwave',nfish);

roikind = 'circle';
% roikind = 'anat';

cond_codes =[400:404];

% PARAMETERS

window.min = [-100 100];
window.max = [0 600];
window.movsum = 50;
window.baseline = [-300 0];

noise.fr_abovenoise = 30;
noise.SDfactor = 4; 

lat_limit = 1000;

method = 'movsum';

%% other settings
Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
Lroi= [2 4 6 8 10 12 14];

labelfilt = 'f309'; %@ SET 
temp = labelfilt(end-2:end);

switch roikind
    case 'circle'
        waves = VSDroiTS.(strcat('circ_filt',temp)).data; % automatically choses the dataset according to the number of filter
    case 'anat'
        waves = VSDroiTS.(strcat('filt',temp)).data; %@ SET ...
end


setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET


pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/boxplot+multcomp/preS_F0_filt09';

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
    latency_out{nroi} = [];
    
    for triali = makeRow(sel_trials)
        
        wave = squeeze(waves(:, nroi,triali));
        output = devo_peak2peak(wave, VSDI.timebase, window,noise, method, 0);
        
        measures.peak2peak(nroi,triali) = output.p2p_value;
        measures.peakminusbasel(nroi,triali) = output.peakminusbasel;
        %                         frames.peaklat(rowi,coli,triali) = output.peaklat_ms;
        %                         frames.p2plat(rowi,coli,triali) = output.p2plat_ms;
        %                         frames.onset30_latency_ms(rowi,coli,triali) = output.onset30_latency_ms;
        measures.onsetnoise_ms(nroi,triali) = output.onsetnoise_ms;
        
        measures.noisethresh(nroi,triali) = output.noisethresh;
        % store trials that will be rejected from latency means
        if measures.onsetnoise_ms(nroi,triali) > lat_limit
            latency_out{nroi} = [latency_out{nroi}  triali];
        end
        %         measures.noisethresh(nroi,triali) = output.noisethresh;
        
        
        
        waveW = wave(output.peakidx(1):output.peakidx(2));
        waveslope = diff(waveW);
        meanslope = mean(waveslope);
        
        measures.meanslope(nroi,triali) = meanslope;
        
        clear output wave waveW waveslope meanslope
        
    end %triali
end %nroi
% 
% subplot(1,5,1)

%% CALCULATE AVERAGE MEASURES 
... the measure has to be computed for each condition label, while for the plots it's calculated trial-wise (and then the boxplot organizes it)

j=1; %
for code = makeRow(cond_codes)
    tricond = intersect(find(VSDI.condition(:,1) == code) , sel_trials);
    
    avewave = mean(waves(:,:,tricond),3) ;

    
    for nroi = makeRow(Rroi)
        roiwave = squeeze(avewave(:,nroi));
        output = devo_peak2peak(roiwave, VSDI.timebase, window, noise, method, 0);
        
        waveW = roiwave(output.peakidx(1):output.peakidx(2));
        slopeval = mean(diff(waveW));
        
        avemeasures.peakminusbasel(nroi,j) = output.peakminusbasel;
        avemeasures.onsetnoise_ms(nroi,j) = output.onsetnoise_ms;
        avemeasures.meanslope(nroi,j) = slopeval;
        
        clear output roiwave waveW slopeval
    end %roi 
    
    % get mA corresponding to each condition (for later plotting)
    mAcond(j)=VSDI.condition(tricond(1),4);
    j = j+1;
end

%% MEASURE SELECTION

% close all
roiname = 'dm4m_R';
% result = 'peakminusbasel';
result ='onsetnoise_ms';
% result = 'meanslope';

boxplot_ylim = [-.01 .04];

% close all

%% BOXPLOT + STAT 
nroi = name2idx(roiname, VSDI.roi.labels);

measure = []; % reinitialize variables
avemeasure = [];

extratitle = '';

switch result
    case 'peakminusbasel'
        measure = makeCol(squeeze(measures.peakminusbasel(nroi,sel_trials)));
        avemeasure = squeeze(avemeasures.peakminusbasel(nroi,:));
        
    case 'meanslope'
        measure = makeCol(squeeze(measures.meanslope(nroi,sel_trials)));
        avemeasure = squeeze(avemeasures.meanslope(nroi,:));
        
    case 'onsetnoise_ms'
        sel_trials = setdiff(sel_trials, latency_out{nroi});
        measure = makeCol(squeeze(measures.onsetnoise_ms(nroi,sel_trials)));
        extratitle = ['-' num2str(length(latency_out{nroi})) 'trials out'];
        
        avemeasure = squeeze(avemeasures.onsetnoise_ms(nroi,:));
        
end

% BUILD MATRIX FOR ANOVA
mA= [];
mA = VSDI.condition(sel_trials,4); %mA

% BOXPLOT
figure
subplot(1,3,1)
boxplot(measure, mA, 'Colors', 'k')
% ylim(boxplot_ylim);

xlabel('mA'); ylabel('% /Delta F');



% subplot(1,3,2) %plot the R-hemisph in the last empty plot (to have a visual guide)
% temp=  squeeze(squeeze(mean(waves(:,nroi,sel_trials),2)));
% plot(VSDI.timebase,temp'); hold on %plots the mean wave of all rois (GS from selected rois)
% %     ylim(local_plotlim)



%STATS

[p, table, stats] = anova1(measure,mA, 'off')

subplot(1,2,1) 
scatter(avemeasure, mAcond, 'filled')
ylabel('mA')
ylim([-0.1 mAcond(end)*1.1])
% xlim([min(avemeasure)*0.9 max(avemeasure)*1.1])
xlim([0 100])
yticks(mAcond)
set(gca, 'ydir', 'reverse')
title('result from averaged wave')



h = subplot(1,2,2)
[c, m , h, gnames] = multcompare_modif(stats, 'CType', 'bonferroni')
title(['multicomp.' extratitle])

sgtitle([num2str(VSDI.ref),'.', VSDI.roi.labels{nroi}, ':', result, '(', labelfilt, ')'])


% end stats

%% LAST UPDATE: 17/05/21
% CREATED: 17/05/21