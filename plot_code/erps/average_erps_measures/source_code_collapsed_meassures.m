%% ERPs

%% CHECK: WAVE MEASURES
for nfish = 1

clearvars -except nfish

reject_on =0;

    % user_settings
    path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/VSDI_ourToolbox/';
    path.TORus = '/home/tamara/Documents/MATLAB/VSDI';
    path.data = fullfile(path.TORus, 'data');
    path.grouplist = path.TORus;
    path.list =fullfile(path.TORus, 'data','BVlists');
    addpath(genpath(path.rootpath));
    addpath(path.TORus);
    % end of user_settings


VSDI = TORus('load',nfish);

% load waves to plot:
VSDroiTS= TORus('loadwave',nfish);

% field = 'circ_filt309' ;
field = 'circ_filt309' ;

waves = VSDroiTS.(field).data; %@ SET: whether circ_filt06, circ_filt09, or filt06

fcode =  [field(1) field(end-2:end)]; fcode = erase(fcode, 'f');

selroi = [3 7];   %@ SET : SELECT ONLY 2 ROIS (the code is hardwired to accept only 2 of them)
% selroi = [7 11];   %@ SET;

setting.manual_reject = 0; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 0; %@ SET

%% SELECT CONDITIONS

% cond_codes =sort(unique(VSDI.condition(:,1)));
% cond_codes = cond_codes(~isnan(cond_codes));
cond_codes = [300:303]; %@SET

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


%% CALCULATE MEASURE FOR THE AVERAGE

window.min = [-100 100];
window.max = [0 600];
window.movsum = 50; %ms
window.slope = 50; %ms

window.baseline = [-300 0];
% window.baseline = [-50 0];

noise.fr_abovenoise = 30;
noise.SDfactor =  4; 

method = 'movsum';

%% ONSET OF  WAVES FROM CIRCULAR ROIS FROM ONE CONDITION AVERAGE
pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/average_erps_measures/intensos';

% build color matrix
colors = lines(length(cond_codes));

figure
ploti = 1;%counter

for nroi = makeRow(selroi)
    
    subplot(2,2,ploti) %waves from all conditions in the first column (loop through conditions)
    hold on;
    
    for codi = 1:length(cond_codes)
        code = cond_codes(codi);
        tricond = intersect(find(VSDI.condition(:,1) == code) , sel_trials);
        codemA(codi)= VSDI.condition(tricond(1),4); %mA corresponding to the code label (for the next subplot)
        roiwave = mean(waves(:,nroi,tricond),3) ;
        
        hold on;
        output = devo_peak2peak(roiwave, VSDI.timebase, window, noise, method, 0);
        
        idx0= dsearchn(VSDI.timebase, 0);
        
        waveW = roiwave(idx0:output.peakidx(2));
        slopemean(nroi, codi) = mean(diff(waveW));
        
        slopemax(nroi,codi) = output.slopemax;
        peak(nroi,codi) = output.peakminusbasel;
        onset_ms(nroi,codi)= output.onsetnoise_ms;
        
        plot(VSDI.timebase,roiwave, 'color', colors(codi,:), 'linewidth', 0.8); hold on %  'linewidth', 1.8
        xline(onset_ms(nroi,codi), 'color', colors(codi,:),'linewidth', 0.8);
        %         xline(0,);
        
%         ylim([-0.2 .3])
        ylim([-0.2 .4]) %for high stimuli

        xlim([-300 600])
        ylabel('%\Delta F (trials ave)');
        
        title([ VSDI.roi.labels{nroi}])
        
        codelegend{codi} = num2str(code);
        
        clear output roiwave waveW slopeval
        
    end %codi
    
%     legend (codelegend); 
    hold off
    ploti = ploti+1;

    % initialize text coord
%     ny = 10; % scaling factor (space between lines)
%     nx = 10; %scaling factor (x space between blocks of text)
%     x = 0;
%     y = length(cond_codes)*ny*6; % 'nºconditions' * 'nº lines per condition' . the number of lines is doubled (x2) to leave a space between lines
x=0; y=0;    
    subplot(2,2,ploti) %waves from all conditions in the first column (loop through conditions)

    xlim([0 length(cond_codes)])
    ylim([0 length(cond_codes)])

    set(gca, 'ydir', 'reverse')
    for codi = 1:length(cond_codes)
        
        text(x, y, {['mA=', num2str(codemA(codi))],...
            ['peak=', num2str(round(peak(nroi,codi),3))],...
            ['slope_m_x=', num2str(round(slopemax(nroi,codi)*100,2))]...
            ['onset=',  num2str(onset_ms(nroi,codi))]...
            }, 'color', colors(codi,:), 'fontsize', 5) 

        
        y = y +1; 
        x = x +1;
    end
    
    ploti = ploti+1;
    
end %  roi
    
        
    if reject_on
        sgtitle([num2str(VSDI.ref),'f:', fcode , '.', num2str(noise.SDfactor), 'sd' '(cl)'])
    else
        sgtitle([num2str(VSDI.ref),'f:', fcode , '.', num2str(noise.SDfactor),'sd'])
    end
    
    
    % Save
    if reject_on
        name = strcat(['collapsed-roi:', num2str(selroi),'_',num2str(VSDI.ref),'filt', fcode, 'block' , num2str(cond_codes(1)), '(clean)']);
    else
        name = strcat(['collapsed-', num2str(selroi),'_',num2str(VSDI.ref),'filt', fcode, 'block' , num2str(cond_codes(1))])
    end
    saveas(gcf,fullfile(pathsave,name),'jpg')
    close
    

end % for nfish