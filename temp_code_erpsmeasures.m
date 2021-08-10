%% ERPs
%% CHECK: WAVE MEASURES
for nfish =11:12

clearvars -except nfish

reject_on = 1;
user_settings;

clearvars -except nfish reject_on

VSDI = TORus('load',nfish);


% load waves to plot:
VSDroiTS= TORus('loadwave',nfish);
waves = VSDroiTS.circ_filt306.data; %@ SET: whether circ_filt06, circ_filt09, or filt06
ppcode = 'c06'; %@ SET according to waves!!!

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

cond_codes =sort(unique(VSDI.condition(:,1)));
cond_codes = cond_codes(~isnan(cond_codes));

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

for condi = makeRow(cond_codes) %to make sure that only the conditions of interest are computed)
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

window.baseline = [-300 0];

noise.fr_abovenoise = 30;
noise.SDfactor = 2;

method = 'movsum';

%% ONSET OF  WAVES FROM CIRCULAR ROIS FROM ONE CONDITION AVERAGE
pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/average_erps_measures';

for code = makeRow(cond_codes)
    tricond = intersect(find(VSDI.condition(:,1) == code) , sel_trials);
    
    avewave = mean(waves(:,:,tricond),3) ;
    figure
    ploti = 1;%counter
    
    for nroi = makeRow(Rroi)
        roiwave = squeeze(avewave(:,nroi));
        output = devo_peak2peak(roiwave, VSDI.timebase, window, noise, method, 0);
        
        waveW = roiwave(output.peakidx(1):output.peakidx(2));
        slopeval = mean(diff(waveW));
        
        subplot(3,3,ploti);
        plot(VSDI.timebase,roiwave); hold on
        xline(output.onsetnoise_ms);
        yline(output.noisethresh, 'r--');
        ylim([-0.2 .3])
        ylabel('%\Delta F (trials ave)');
        text(output.onsetnoise_ms,-0.1,num2str(output.onsetnoise_ms),'FontSize',6)
        
        text(-500,0.2,['A=' num2str(round(output.p2p_value,2))],'FontSize',6) % peak-to-peak amplitude
        text(-500,0.15,['S=' num2str(round(slopeval*100,2))],'FontSize',6) %

        title([ VSDI.roi.labels{nroi}])
        ploti = ploti+1;
        
        clear output roiwave waveW slopeval
    end %roi 
   
    ...mute if the roi are not circular
    ax9= subplot(3,3,9);
    roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(Rroi,:), VSDI.roi.circle.R, ax9);
    axis image
    
    if reject_on
        sgtitle([num2str(VSDI.ref),'pp', ppcode , 'average waves cond=', num2str(code), '(cl)'])
    else
        sgtitle([num2str(VSDI.ref),'pp', ppcode , 'average waves cond=', num2str(code)])
    end
    
    
    % Save
    if reject_on
        name = strcat([num2str(VSDI.ref),'pp', ppcode , 'average waves cond=', num2str(code),'(clean)']);
    else
        name = strcat([num2str(VSDI.ref),'pp', ppcode , 'average waves cond=', num2str(code)])
    end
    saveas(gcf,fullfile(pathsave,name),'jpg')
    close
    
end

end ... for nfish