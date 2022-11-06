%% MEAN BRADYCHARDIA GROUPED AND EXCEL OUTPUT WITH AMPLITUDE MEASURES
% The bradychardia is calculated leaving out a window with the US artifact
% (when it's tone, indicated by VSDI.condicion(:,1) == NaN, no safety window is applied)

% Ampitude measures from main rois are also reflected

% Shark/noShark classification is made such as if dldm's peak is higher than
% dm4'speak, a shark is considered. That peak is to be found in the whole trial, unlike the peaks used to measure the rois peaks
...that are to be found according to the input 'window'
    
clear

%SETTINGS
export.excel = 1; % '1' = sí; '0'=no
w.pre= 10; % ventana pre-estímulo(en segundos) , de los que se va a tomar las espigas para hacer la media de frecuencia instantánea o el conteo
w.post = w.pre/2; %sólo para el conteo
outfolder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady';
excelname = fullfile(outfolder,['brady_grouped_bycondition.xlsx']);

% set(0,'DefaultFigureVisible','off')


% SETTINGS FOR THE FUNCTION
window.min = [-100 100];
window.max = [0 1000];
window.movsum = 50;
window.basel = [-100 0];
window.slope=50;


method = 'movsum';
% END OF SETTINGS

% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));
%END OF USER SETTINGS

%%
for nfish = [11 12]
    
    VSDI = TORus('load',nfish);
        spike = TORus('loadspike', nfish); % ECG
    VSDroiTS =TORus('loadwave',nfish);
    waves = VSDroiTS.circ_filt309.data; %@ SET
    
    
    roi2plotidx = [1 3 5 7 11];  %@ SET
    
    window_sh= window; % DIFFERENT WINDOW FOR SHARK-FINDING  (the whole trial). It has to be defined inside the function because VSDI.timebase might differ among fishes
    window_sh.max = [0 VSDI.timebase(end-10)];
    
    
    %     output = {'trial', 'spiketime', 'cond','mA','', '%brady(count)','', 'IBIbasel', 'ibi0', 'ibi1','ibi2','ibi3','', 'estimation','dm4','dm4m','dm2','dldm'}; %header first
    
    
    for triali =  makeRow (VSDI.nonanidx)
        
        % ----------------------------------------------
        % SPIKE(ecg)-BASED MEASURES
        % ----------------------------------------------
        if  VSDI.spike.RO(triali) == 0
            flag_noROspike = 1;
        else
            flag_noROspike = 0;
            
        end
        
        
        if ~flag_noROspike
            
            stim = VSDI.spike.RO(triali)+VSDI.info.Sonset/1000 ; %stimuus arrival time (in spike) for this trial (RO time here has already left out the shutter delay, so we dont have to sum it up)
            
            
            % 1. LEAVE OUT STIMULUS NOISE (if it's a somatic stimulus):
            
            % Get safety window to remove US noise
            if (VSDI.list(triali).Sdur >= 60) & (VSDI.condition(triali,4) >= 1) % EXPERIMENT-SPECIFIC PARAMETER: if it is a long train (60ms) of 1mA or more
                noisewin = 430/1000; % (430ms) empirically measured from the noise captured by the cup-electrodes
            else
                noisewin = VSDI.info.Sdur(triali)/1000+ 0.03; %+0.03security window (from spike, the events normally have a 0.02 minimum step from spike to spike)
            end
            
            % 2. Clean the spikes due to that trial's noise
            clean_spikes = spike.spikes(spike.spikes<stim -0.01); %keep all before the stim (excluding another prestim safety window)
            
            if ~isnan(VSDI.condition(triali,1))% IF IT IS NaN, IT IS A TONE, so the artifact due to the stim does not exist
                clean_post = spike.spikes(spike.spikes> stim + noisewin );
            else % If it is a tone, no 'noisewin' clearing is needed
                clean_post = spike.spikes(spike.spikes> stim);
                
            end
            clean_spikes = [clean_spikes; clean_post];
            
            % 3. GET IDX FOR BRADYCHARDIA MEASURE FROM THE CLEANED HEART
            % SPIKES
            idx_pre = find(clean_spikes > stim-w.pre  & clean_spikes < stim);
            
            idx_post =  find(clean_spikes > stim  & clean_spikes < stim+ w.post);
            idx_spike0 = find(clean_spikes> stim, 1, 'first');
            
            % bradycardia counting
            precount = length(idx_pre);
            postcount = length(idx_post)*2;
            brady_count = ((precount - postcount)/precount)*100;
            
            % IBI calculation
            preibi = mean(diff(clean_spikes(idx_pre)));
            
            post0 = clean_spikes(idx_spike0) - clean_spikes(idx_spike0-1);
            post1 = clean_spikes(idx_spike0+1) - clean_spikes(idx_spike0);
            post2 = clean_spikes(idx_spike0+2) - clean_spikes(idx_spike0+1);
            post3 = clean_spikes(idx_spike0+3) - clean_spikes(idx_spike0+2);
            
            % PERCENT IBI calculation
            perc0 = ((post0 - preibi) / preibi)*100;
            perc1 = ((post1- preibi) / preibi)*100;
            perc2 = ((post2 - preibi) / preibi)*100;
            perc3 = ((post3 - preibi) / preibi)*100;
            
        end % flag
        
        % ----------------------------------------------
        % BRAINVISION/WAVE-BASED MEASURES
        % ----------------------------------------------
        for nroi = roi2plotidx
            wave = squeeze(waves(:, nroi,triali));
            local_output = devo_peak2peak(wave, VSDI.timebase, window,[], method, 0);
            local_output_sh = devo_peak2peak(wave, VSDI.timebase, window_sh,[], method, 0);
            
            peak(nroi) = local_output.peakminusbasel;
            peak_sh(nroi) = local_output_sh.peakminusbasel;
        end
        
        
        % ----------------------------------------------
        % OUTPUT FOR EXCEL
        % ----------------------------------------------
        
        
        if ~flag_noROspike
            output(triali,6) = round(brady_count);
            
            output(triali,8) = round(preibi,2);
            output(triali,9) = round(post0,2);
            output(triali,10) = round(post1,2);
            output(triali,11) = round(post2,2);
            output(triali,12) = round(post3,2);
            
            output(triali,14) = round(perc0,2);
            output(triali,15) = round(perc1,2);
            output(triali,16) = round(perc2,2);
            output(triali,17) = round(perc3,2);
            
            
        end
        
        output(triali,20) = round(peak(1),3); %dm4
        output(triali,21) = round(peak(3),3); %dm4m
        output(triali,22) = round(peak(7),3); % dm2
        output(triali,23) = round(peak(11),3);% dldm
        
        
        clearvars -except VSDI spike VSDroiTS waves roi2plotidx output j export w outfolder window method window_sh nfish excelname
        
    end %triali
    
    % ----------------------------------------------
    % MEAN VALUES FOR EACH CONDITION
    % ----------------------------------------------
    trial_kinds = unique(VSDI.condition(:,1));
    trial_kinds = trial_kinds(~isnan(trial_kinds));
    
    % HEADER
    localoutput{1,1} = 'cod';
    localoutput{1,2} = 'mA';
    
    localoutput{1,3} = '%brady (beats)'; %brady count
    
    localoutput{1,5} = '%ibi0'; % ibis0 %
    localoutput{1,6} = '%ibi1'; % ibis1 %
    localoutput{1,7} = '%ibi2'; % ibis2 %
    localoutput{1,8} = '%ibi3'; % ibis3 %
    
    localoutput{1,10} ='dm4'; % dm4m
    localoutput{1,11} = 'dm4m';  % dm4m
    localoutput{1,12} ='dm2';  % dm2
    localoutput{1,13} = 'dldm'; % dldm
    
    
    k = 2;
    for which_cond =  makeRow(trial_kinds)
        
        localoutput{k,1} = which_cond;
        
        
        cond_trials = find(VSDI.condition(:,1) == which_cond);
        
        localoutput{k,2} = VSDI.condition(cond_trials(1),4);
        
        localoutput{k,3} = round(mean(output(cond_trials,6)),2); %brady count
        
        localoutput{k,5} = round(mean(output(cond_trials,14)),2); % ibis0 %
        localoutput{k,6} = round(mean(output(cond_trials,15)),2); % ibis1 %
        localoutput{k,7} = round(mean(output(cond_trials,16)),2); % ibis2 %
        localoutput{k,8} = round(mean(output(cond_trials,17)),2); % ibis3 %
        
        localoutput{k,10} = round(mean(output(cond_trials,20)),3); % dm4
        localoutput{k,11} = round(mean(output(cond_trials,21)),3); % dm4m
        localoutput{k,12} = round(mean(output(cond_trials,22)),3); % dm2
        localoutput{k,13} = round(mean(output(cond_trials,23)),3); % dldm
        
%         % for the plot
%         condlabel(k-1) = which_cond;
%         roiactiv(k-1,1) = mean(output(cond_trials,20)); % dm4
%         roiactiv(k-1,2) = mean(output(cond_trials,21)); % dm4m
%         roiactiv(k-1,3) = mean(output(cond_trials,22)); % dm2
%         roiactiv(k-1,4) = mean(output(cond_trials,23)); % dldm

        
        k = k+1;

    end
    
%     figure
%     plot(1:length(condlabel),roiactiv(:, 1)); 
%     hold on
%     plot(roiactiv(:, 2)); 
%     plot(roiactiv(:, 3)); 
%     plot(roiactiv(:, 4));
%     legend dm4 dm4m dm2 dldm
%     
%     xticklabels(condlabel)
%     
%     title ([num2str(VSDI.ref), ': peak minus baseline'])

%     ----------------------------------------------
%     WRITE EXCEL
%     ----------------------------------------------
    if export.excel == 1
        % write output (new sheet for each fish
        localoutput = cell2table(localoutput);
        writetable (localoutput, excelname, 'sheet', num2str(VSDI.ref))
    end
    
    
    clearvars -except  export w outfolder window method nfish excelname
    
end %nfish

%% WRITE PARAMETERS USED IN ANALYSIS

if export.excel == 1
    % write output (new sheet for each fish
    paramet{1,1} = 'cod';
    
    paramet{1,1}= 'mA ';
    paramet{1,2}= w.pre;
    paramet{1,3}= 'postStim window (s)';
    paramet{1,4}= w.post;
    
    paramet{1,5}= 'Peak minus baseline measure is shown for each roi';
    
    paramet{6,1}= 'window for peak finding:'; paramet{6,2}= [num2str( window.max(1)) 'to' num2str( window.max(2))  'ms'];
    paramet{7,1}= 'baseline:'; paramet{7,2}=  [num2str(window.basel(1)) 'to' num2str(window.basel(2))  'ms'];
    
    paramet{9,1} = 'time window for shark findind spans the whole trial (unlike for roi peak activity)';
    
    
    paramet =  cell2table(paramet);
    writetable (paramet, excelname, 'sheet', 'param')
end

blob()



%% PLOT RELATIONSHIP OF BRADYCHARDIA AND CONDITIONS IN EACH


%% Created: 14/07/2021

