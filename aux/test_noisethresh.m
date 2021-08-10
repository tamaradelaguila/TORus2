
clear
user_settings
nfish =  12%@ SET

[VSDI] = TORus('load',nfish);

VSDroiTS= TORus('loadwave',nfish);
waves = VSDroiTS.circ_filt306.data; %@ SET


% temp = TORus('loadmovie',nfish,'_06filt3');
% movies = temp.data(:,:,1:end-1,:);
idx = VSDI.nonanidx(1:16);
% idx = setdiff(idx, VSDI.reject.GSabs025); 
% idx = setdiff(idx, VSDI.reject.GSdeviat2sd); 


% function settings

window.min = [-100 100];
window.max = [0 600];
window.movsum = 50;

window.baseline = [-100 0];
baseidx = find_closest_timeidx(window.baseline, VSDI.timebase);

noise.fr_abovenoise = 30;
noise.SDfactor = 3; 

ylimite = [-.05 0.20];

nroi= 3; 
%%
figure
sgtitle(['noise threshold + onset (>mean +' num2str(noise.SDfactor) 'Std in' num2str(noise.fr_abovenoise) 'consec.frames). Basel='  num2str(window.baseline(1)) 'to'  num2str(window.baseline(2)) 'ms'])
for ploti = 1:15
    subplot(4,4,ploti)
    
%     wave = squeeze(movies(15,15,:,idx(ploti))); % from movie
wave = squeeze(waves(:,nroi,idx(ploti))); 
    
    output = devo_peak2peak(wave, VSDI.timebase, window, noise, 'movsum', 0);
    
    plot(VSDI.timebase, wave, 'linewidth', 1.5); hold on;
    yline(output.noisethresh, 'g', 'linewidth', 1.5)
    xline(output.onsetnoise_ms, 'r--', 'linewidth', 1.5)
    ylim(ylimite)
    xline(window.baseline(1)); xline(window.baseline(2));   
    yline(mean(wave(baseidx(1):baseidx(2))), 'b');
%     legend data noise_thresh onset
    clear wave  output
end

ploti =16
    subplot(4,4,ploti)
    
    wave = squeeze(waves(:,nroi,idx(ploti))); 

%     wave = squeeze(movies(15,15,:,idx(ploti)));
    
    output = devo_peak2peak(wave, VSDI.timebase, window, noise, 'movsum', 0);
    
    plot(VSDI.timebase, wave, 'linewidth', 1.5); hold on;
    yline(output.noisethresh, 'g', 'linewidth', 1.5)
    xline(output.onsetnoise_ms, 'r--', 'linewidth', 1.5)
    ylim(ylimite)
    xline(window.baseline(1)); xline(window.baseline(2));   

    legend data noise_thresh onset
    clear wave  output

   %% PLOT BY CONDITION
   cond_code =[402 ];
setting.manual_reject =1;
% SELECT EXCLUDED

rejectidx = [];

    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];


    rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];
    


    rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];
    rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
    

rejectidx = sort(unique(rejectidx));
    

% SELECT CASES 
    
        sel_trials = makeCol(find(VSDI.condition(:,1)==cond_code));
    
        sel_trials = setdiff(sel_trials, rejectidx);
    
    
%%
figure
sgtitle(['noise threshold + onset (>mean +' num2str(noise.SDfactor) 'Std in' num2str(noise.fr_abovenoise) 'consec.frames). Basel='  num2str(window.baseline(1)) 'to'  num2str(window.baseline(2)) 'ms'])
for ploti = 1:15
    subplot(4,4,ploti)
    
%     wave = squeeze(movies(15,15,:,idx(ploti))); % from movie
wave = squeeze(waves(:,nroi,sel_trials(ploti))); 
    
    output = devo_peak2peak(wave, VSDI.timebase, window, noise, 'movsum', 0);
    
    plot(VSDI.timebase, wave, 'linewidth', 1.5); hold on;
    yline(output.noisethresh, 'g', 'linewidth', 1.5);
    xline(output.onsetnoise_ms, 'r--', 'linewidth', 1.5)
    ylim(ylimite)
    xline(window.baseline(1)); xline(window.baseline(2));   
    yline(mean(wave(baseidx(1):baseidx(2))), 'b');
%     legend data noise_thresh onset
    clear wave  output
end

ploti =16
    subplot(4,4,ploti)
    wave = squeeze(waves(:,nroi,sel_trials)); 
    erps = mean(squeeze(waves(:,2,sel_trials)),2); 

%     wave = squeeze(movies(15,15,:,idx(ploti)));
    
    output = devo_peak2peak(erps, VSDI.timebase, window, noise, 'movsum', 0);
    
    % PLOT
    plot(VSDI.timebase, wave); hold on;

    plot(VSDI.timebase, erps, 'k', 'linewidth', 1.5); hold on;
    yline(output.noisethresh, 'g', 'linewidth', 1.5); 
    xline(output.onsetnoise_ms, 'k--', 'linewidth', 1.5)
    ylim(ylimite)
    xline(window.baseline(1)); xline(window.baseline(2));   
    title('erps')
    display(output.onsetnoise_ms)
    legend data noise_thresh onset
    clear wave  output
    
    
    
    
    
%% PLOT AND PRINT  LATENCY RESULTS
    
% ----------------------------------------------
% SETTINGS
% ----------------------------------------------

clear
close all

nfish = 12;
reject_on = 1;
user_settings;

%load structure + wave
VSDI = TORus('load',nfish);
VSDroiTS= TORus('loadwave',nfish);

roikind = 'circle';

cond_codes =[400:404];

% PARAMETERS

window.min = [-100 100];
window.max = [0 600];
window.movsum = 50;
window.baseline = [-300 0];

noise.fr_abovenoise = 30;
noise.SDfactor = 3; 

method = 'movsum';


lat_limit = 1000;
% ----------------------------------------------
% SELECT DATA 
% ----------------------------------------------
switch roikind
    case 'circle'
        waves = VSDroiTS.circ_filt306.data; %@ SET
    case 'anat'
        waves = VSDroiTS.filt306.data; %@ SET
end



Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
% Lroi= [2 4 6 8 10 12 14];

% ----------------------------------------------
% SELECT EXCLUDED
% ----------------------------------------------

setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET

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

% ----------------------------------------------
% SELECT CASES
% ----------------------------------------------

sel_trials= [];

for condi = cond_codes %to make sure that only the conditions of interest are computed)
    condtrials = makeCol(find(VSDI.condition(:,1)==condi));
    sel_trials  = [sel_trials; condtrials];
end
sel_trials= sort(sel_trials);

if reject_on  %@ SET
    sel_trials = setdiff(sel_trials, rejectidx);
end

% ----------------------------------------------
% CALCULATE MEASURE FOR EACH TRIAL
% ----------------------------------------------

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


% --------------------------------------------------------
% PLOT & SAVE -FROM EACH TRIAL- ALL ROIS FROM ONE HEMISPHERE FOR EACH TRIAL
% --------------------------------------------------------

path_latencyplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_noisethresh/test_threshold'; %@ SET

% FIND YLIM VALUES AMONG NON-OUTLIER MAX AND MIN VALUES
tempwave = squeeze(waves(:, nroi,sel_trials));
tempmax = max(tempwave); 
tempmax = tempmax(~isoutlier(tempmax)); %select only maximum values that are not outliers
tempmax = max(tempmax);

tempmin = min(tempwave); 
tempmin = tempmin(~isoutlier(tempmin)); %select only maximum values that are not outliers
tempmin = mean(tempmin);


ylimit = [tempmin tempmax];

clear  tempwave tempmax tempmin

for triali =  makeRow(sel_trials)
    
    ploti = 1;
    for nroi = Rroi(1:end-1)
        subplot(3,3,ploti)
        wave = squeeze(waves(:, nroi,triali));
        thresh = measures.noisethresh(nroi,triali);
        latency= measures.onsetnoise_ms(nroi,triali);
        plot(VSDI.timebase, wave);
        yline(thresh, 'b'); xline(latency, 'g--', 'linewidth',1.2);
        title (VSDI.roi.labels{nroi});
        ploti = ploti+1;
        ylim (ylimit);
    end
    
    subplot(3,3,7) %last plot apart to draw legend
    nroi = Rroi(end);
    wave = squeeze(waves(:, nroi,triali));
    thresh = measures.noisethresh(nroi,triali);
    latency= measures.onsetnoise_ms(nroi,triali);
    
    plot(VSDI.timebase, wave);
    yline(thresh, 'b'); xline(latency, 'g--', 'linewidth',1.2);
    title  (VSDI.roi.labels{nroi});
    ylim (ylimit);
    
    legend data noise-thresh latency
    
    axH = subplot(3,3,9);
    roicirc_preview_multiple(VSDI.backgr(:,:,triali), VSDI.roi.circle.center, VSDI.roi.circle.R, axH)
    
    trialname = num2str(VSDI.trialref(triali)); trialname = [trialname(end-2:end) 'A'];
    name = [num2str(VSDI.ref), '.trial=', trialname, '-(', num2str(VSDI.condition(triali,4)), 'mA)'];

        sgtitle(name)
    
    
    
    saveas(gca, fullfile(path_latencyplot, [name '.jpg']), 'jpg')
    close
    clear wave name trialname
end

% --------------------------------------------------------
% PLOT & SAVE ALL ROIS FROM ONE HEMISPHERE FOR THE ERPS
% --------------------------------------------------------
ylimit = [-0.05 0.25]; 
for cond_codes =[400:404]
    cond_trials = find(VSDI.condition(:,1) == cond_codes);
    cond_trials =  intersect(cond_trials, sel_trials);
    
    ploti = 1;
    for nroi =  Rroi(1:end-1)
        out = latency_out{nroi};
        subplot(3,3,ploti)
        latentrial = setdiff(cond_trials, out); % use data about rejected trials  for the current roi (stored from the trial-wise result of devo_peak2peak)
%         latentrial = cond_trials;
        wave = mean(waves(:, nroi,latentrial),3);
        
        erps = devo_peak2peak(wave, VSDI.timebase, window,noise, method, 0); %calculate the measures for the erp
        
        plot(VSDI.timebase,wave);
        yline(erps.noisethresh, 'b'); xline(erps.onsetnoise_ms, 'g--', 'linewidth',1.2);
        
        ylim (ylimit);
        
        title = (VSDI.roi.labels{nroi});
        ploti = ploti+1;
        clear wave erps thresh latency
    end
    
    subplot(3,3,7) %last plot apart to draw legend
    nroi = Rroi(end);
    latentrial = setdiff(cond_trials, out); % use data about rejected trials  for the current roi (stored from the trial-wise result of devo_peak2peak)
    
    wave = mean(waves(:, nroi,latentrial),3);
    
    erps = devo_peak2peak(wave, VSDI.timebase, window,noise, method, 0); %calculate the measures for the erp
    
    plot(VSDI.timebase,wave);
    yline(erps.noisethresh, 'b'); xline(erps.onsetnoise_ms, 'g--', 'linewidth',1.2);
    ylim (ylimit);
    
    title = (VSDI.roi.labels{nroi});
    
    legend data noise-thresh latency
    
    axH = subplot(3,3,9);
    roicirc_preview_multiple(VSDI.backgr(:,:,triali), VSDI.roi.circle.center, VSDI.roi.circle.R, axH)
    
    trialname = num2str(VSDI.trialref(triali)); trialname = [trialname(end-2:end) 'A'];
    name = [num2str(VSDI.ref),'-', '.erps from condition=', num2str(cond_codes)];
    sgtitle(name)
    
    saveas(gca, fullfile(path_latencyplot, [name '.jpg']), 'jpg')
    close
    clear wave name  erps
end

