%% FIGURE SKETCH - 

%% C: mA vs peak


clear
% clearvars -except condition_list

%SETTINGS
outfolder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady';
excelname = fullfile(outfolder,'brady_grouped_bycondition.xlsx');

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

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/02_seleccion_indiv_agrupadosBradi/plots_amplitude' ;%@ SET

%END OF USER SETTINGS
roi2plotidx = [1 3 5 7 11];  %@ SET
load('/home/tamara/Documents/MATLAB/VSDI/TORus/data/condition_list.mat')

% ---------------------------------------------------------------------------
%% AMPLITUDE TRIAL-WISE TO PLOT BY CONDITION
% ---------------------------------------------------------------------------

for  block = 16:17%1:length(condition_list)
    
    nfish = condition_list{block,1};
    
    trial_kinds = condition_list{block,2};
    cond_def =condition_list{block,3};
    
    VSDI = TORus('load',nfish);
    %     spike = TORus('loadspike', nfish); % ECG
    VSDroiTS =TORus('loadwave',nfish);
    waves = VSDroiTS.circ_filt309.data; %@ SET
    
    
    %     output = {'trial', 'spiketime', 'cond','mA','', '%brady(count)','', 'IBIbasel', 'ibi0', 'ibi1','ibi2','ibi3','', 'estimation','dm4','dm4m','dm2','dldm'}; %header first
    
    
    for triali =  makeRow (VSDI.nonanidx)
        
        % ----------------------------------------------
        % BRAINVISION/WAVE-BASED MEASURES
        % ----------------------------------------------
        for nroi = roi2plotidx
            wave = squeeze(waves(:, nroi,triali));
            local_output = devo_peak2peak(wave, VSDI.timebase, window,[], method, 0);
            
            peak(triali, nroi) = local_output.peakminusbasel;
            clear wave local_output
        end
        
        
        
    end %triali
    
    % ----------------------------------------------
    % MEAN VALUES FOR EACH CONDITION
    % ----------------------------------------------
    
    % HEADER
    j = 1;
    for condi = makeRow( trial_kinds)
        sel_trials = find(VSDI.condition(:,1) == condi);
        meanAct(j,1) = mean(peak(sel_trials,1)); %dm4
        meanAct(j,2) = mean(peak(sel_trials,3)); %dm4m
        meanAct(j,3) = mean(peak(sel_trials,7)); %dm2
        meanAct(j,4) = mean(peak(sel_trials,11)); %dldm
        
        stdAct(j,1) = std(peak(sel_trials,1)); %dm4
        stdAct(j,2) = std(peak(sel_trials,3)); %dm4m
        stdAct(j,3) = std(peak(sel_trials,7)); %dm2
        stdAct(j,4) = std(peak(sel_trials,11)); %dldm
        
        mA(j) = VSDI.condition(sel_trials(1),4) ;
        
        j=j+1;
        clear sel_trials
    end
    
    %         % for the plot
    %         condlabel(k-1) = which_cond;
    %         roiactiv(k-1,1) = mean(output(cond_trials,1)); % dm4
    %         roiactiv(k-1,2) = mean(output(cond_trials,2)); % dm4m
    %         roiactiv(k-1,3) = mean(output(cond_trials,3)); % dm2
    %         roiactiv(k-1,4) = mean(output(cond_trials,4)); % dldm
    
    if mA(2) == 0 || isnan(mA(2)) % if the second condition is NaN or 0, then it is a tone and some adjustements are needed to be plotted
        flag_adapt =1;
        mA(1) = -0.5; mA(2) = 0;
        
        labels{1} = '0';labels{2} = 'tone';
        
        for ii = 3:length(mA)
            labels{ii} =  num2str(mA(ii));
        end
    else
        
        flag_adapt =0;
        
    end
    
    
    figure
    errorbar(mA, meanAct(:,1), stdAct(:,1), 'o-','linewidth', 1.8)
    hold on
    errorbar(mA, meanAct(:,2), stdAct(:,2), 'o-','linewidth', 1.8)
    errorbar(mA, meanAct(:,3), stdAct(:,3), 'o-','linewidth', 1.8)
    errorbar(mA, meanAct(:,4), stdAct(:,4), 'o-','linewidth', 1.8)
    legend dm4 dm4m dm2 dldm
    
    xticks(mA);
    
    if flag_adapt
        xticklabels(labels)
    end
    
    
    xlim([min(mA)-0.5, max(mA)+0.5])
    xlabel('mA'); ylabel('peak minus baseline')
    title (['peak-B:' num2str(VSDI.ref), cond_def, '(', num2str(trial_kinds(1)) , ')'])
    
    % ----------------------------------------------
    % WRITE EXCEL
    % ----------------------------------------------
    %     if export.excel == 1
    %         % write output (new sheet for each fish
    %         localoutput = cell2table(localoutput);
    %         writetable (localoutput, excelname, 'sheet', num2str(VSDI.ref))
    %     end
    
    trialref = num2str(VSDI.trialref(triali)); trialref = [trialref(end-2:end) 'A'];
    name2save = fullfile(savein,['plot_Amplitude_trialwise' num2str(VSDI.ref) num2str(condi(1)) '.jpg']);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
    saveas(gcf,name2save,'jpg')
    
    close
    clear mA meanAct stdAct  peak
    
end %nfish

blob()