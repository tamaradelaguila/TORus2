%% CHECK: WAVE MEASURES
clear 

reject_on = 1;
user_settings;
nfish =12;

clearvars -except nfish reject_on

VSDI = TORus('load',nfish);


% load waves to plot:
VSDroiTS= TORus('loadwave',nfish);

% waves = VSDroiTS.circ_filt306.data; %@ SET
waves = VSDroiTS.circ_filt309.data; %@ SET...
label = 'c09';%@ ...and SET accordingly

Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
Lroi= [2 4 6 8 10 12 14];

setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET


%% SELECT CONDITIONS
% cond_codes = unique(VSDI.condition(:,1));
% cond_codes=  cond_codes(~isnan(cond_codes));
% cond_codes= setdiff(cond_codes,0);

cond_codes =[400:404 ];

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
        
        method = 'movsum';
    
        
%% WAVES FROM CIRCULAR ROIS FROM ONE CONDITION AVERAGE
noise.SDfactor = 4; 
cond_code= 404;
tricond = intersect(find(VSDI.condition(:,1) == cond_code) , sel_trials);

avewave = mean(waves(:,:,tricond),3) ;
figure
ploti = 1;%counter
for nroi = makeRow(Rroi)
    roiwave = squeeze(avewave(:,nroi));
    output = devo_peak2peak(roiwave, VSDI.timebase, window,[], method, 0);
    
    subplot(3,3,ploti);
    plot(VSDI.timebase,roiwave); hold on
    xline(output.onsetnoise_ms);
    yline(output.noisethresh, 'r--');
    ylim([-0.2 .4])
    ylabel('%\Delta F (trials ave)');
    title([ VSDI.roi.labels{nroi}])
    ploti = ploti+1;
    
    clear output roiwave
end %roi

    ax9= subplot(3,3,9); 
    roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(Rroi,:), VSDI.roi.circle.R, ax9); 
    axis image
    
    sgtitle([num2str(VSDI.ref), 'average waves cond=', num2str(cond_code)])
    
%% WAVES FROM CIRCULAR ROIS FROM ONE SINGLE TRIAL
noise.SDfactor = 2; 
     
triali = 102;
wave = squeeze(waves(:,:,triali)) ;
figure
ploti = 1;%counter
for nroi = makeRow(Rroi)
    roiwave = squeeze(wave(:,nroi));
    output = devo_peak2peak(roiwave, VSDI.timebase, window,noise, method, 0);
    
    subplot(3,3,ploti);
    plot(VSDI.timebase,roiwave); hold on
    xline(output.onsetnoise_ms);
    yline(output.noisethresh, 'r--');
    ylim([-0.2 .5])
    ylabel('%\Delta F (trials ave)');
    title([ VSDI.roi.labels{nroi}])
    ploti = ploti+1;
    
    clear output roiwave
end %roi

    ax9= subplot(3,3,9); 
    roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(Rroi,:), VSDI.roi.circle.R, ax9); 
    axis image
    tempref = VSDI.list(triali).Name ; tempref = tempref(end-7:end-4);
    sgtitle([num2str(VSDI.ref), ' waves trial=',tempref,'(',label,')' ])
