%% FIGURE FAST SKETCH -

clear

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')

% COMMON SETTINGS FOR THE FUNCTION
window.min = [-100 100];
window.max = [0 1000];
window.movsum = 50;
window.basel = [-100 0];
window.slope=50;

noise.fr_abovenoise = 30;
noise.SDfactor = 2;

% Window For average-based analysis

window_ave = window;

noise_ave = noise;
noise_ave.SDfactor = 4;% SET differences


method = 'movsum';
% END OF SETTINGS

% user_settings


% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

% end of user_settings

%% A: TILES CON CURVITA
ref_movie= '_12filt5' ;


%% C: PLOT mA vs peak

clear
% clearvars -except fast_condition_list

%SETTINGS
outfolder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady';
excelname = fullfile(outfolder,'brady_grouped_bycondition.xlsx');

% set(0,'DefaultFigureVisible','off')


%END OF USER SETTINGS
roi2plotidx = [1 3 5 7 11];  %@ SET

% ---------------------------------------------------------------------------
% AMPLITUDE TRIAL-WISE TO PLOT BY CONDITION
% ---------------------------------------------------------------------------

for  block = 1:length(fast_condition_list)
    
    nfish = fast_condition_list{block,1};
    
    trial_kinds = fast_condition_list{block,2};
    cond_def =fast_condition_list{block,3};
    
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
        
        %         stdAct(j,1) = std(peak(sel_trials,1)); %dm4
        %         stdAct(j,2) = std(peak(sel_trials,3)); %dm4m
        %         stdAct(j,3) = std(peak(sel_trials,7)); %dm2
        %         stdAct(j,4) = std(peak(sel_trials,11)); %dldm
        
        semAct(j,1) =  std(peak(sel_trials,1))/sqrt(length(peak(sel_trials,1)));%dm4
        semAct(j,2) =  std(peak(sel_trials,1))/sqrt(length(peak(sel_trials,3))); %dm4m
        semAct(j,3) =  std(peak(sel_trials,1))/sqrt(length(peak(sel_trials,7)));%dm2
        semAct(j,4) =  std(peak(sel_trials,1))/sqrt(length(peak(sel_trials,11)));%dldm
        
        
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
    errorbar(mA, meanAct(:,2), semAct(:,1), 'o-','linewidth', 1.8)
    hold on
    errorbar(mA, meanAct(:,3), semAct(:,2), 'o-','linewidth', 1.8)
    legend dm4m dm2
    
    %     errorbar(mA, meanAct(:,1), semAct(:,3), 'o-','linewidth', 1.8)
    %     errorbar(mA, meanAct(:,4), semAct(:,4), 'o-','linewidth', 1.8)
    %     legend dm4 dm4m dm2 dldm
    
    xticks(mA);
    
    if flag_adapt
        xticklabels(labels)
    end
    
    
    xlim([min(mA)-0.5, max(mA)+0.5])
    ylim([0 0.6])
    xlabel('mA'); ylabel('peak minus baseline mean \pm sem')
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

%% D – CURVAS COND SOLAPADAS POR REGIÓN
%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
for  block = 1:length(fast_condition_list)
    
    nfish = fast_condition_list{block,1};
    cond_codes = fast_condition_list{block,2};
    descript = fast_condition_list{block,3};
    
    %----------------------------------------------------------------
    % @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
    %----------------------------------------------------------------
    measure = {'peakminusbasel'};
    
    %----------------------------------------------------------------
    % @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
    %----------------------------------------------------------------
    reject_on= 0;
    
    roiname1 = 'dm4m_R';
    roiname2 = 'dm2_R';
    
    wavesylim.ave = [-0.2 0.4]; %ylim for pltting average waves
    wavesylim.all = [-0.2 0.6]; %ylim for plotting all waves (has to be higher)
    
    
    [VSDI] = TORus('load',nfish);
    
    tempmov = TORus('loadmovie',nfish,'_09filt3');
    movies = tempmov.data(:,:,1:end-1,:);
    
    display([num2str(VSDI.ref) '-block:' num2str(cond_codes(1)) descript]) % DISPLAY FISH AND BLOCK INFO
    
    VSDroiTS = TORus('loadwave',nfish);
    fcode =  'filt309';
    field = ['circ_' fcode] ;
    waves = VSDroiTS.(field).data; %@ SET: whether circ_filt06, circ_filt09, or filt06
    
    %----------------------------------------------------------------
    % @SET: PEAK2PEAK CONFIGURATION
    %----------------------------------------------------------------
    
    % roiname to roi-idx
    nroi1 = name2idx(roiname1, VSDI.roi.labels);
    nroi2 = name2idx(roiname2, VSDI.roi.labels);
    
    selroi = [nroi1 nroi2];
    
    
    clearvars -except measure reject_on nfish cond_codes VSDI movies waves fcode ...
        window window_ave noise noise_ave rejectidx roiname1 roiname2 method lat_limit...
        wavesylim selroi setting savein fast_condition_list descript
    
    %@ SET : SELECT ONLY 2 ROIS (the code is hardwired to accept only 2 of them)
    
    % SELECT CASES
    sel_trials= [];
    
    for condi = makeRow(cond_codes) %to make sure that only the conditions of interest are computed)
        condtrials = makeCol(find(VSDI.condition(:,1)==condi));
        sel_trials  = [sel_trials; condtrials];
    end
    
    sel_trials= sort(sel_trials);
    
    
    % ONSET OF  WAVES FROM CIRCULAR ROIS FROM ONE CONDITION AVERAGE
    
    % build color matrix
    colors = lines(length(cond_codes));
    
    figure
    ploti = 1;%counter
    
    for nroi = makeRow(selroi)
        
        subplot(1,2,ploti) %waves from all conditions in the first column (loop through conditions)
        hold on;
        
        j = 1; %counter for legend labels
        
        for codi = 1:length(cond_codes)
            code = cond_codes(codi);
            tricond = intersect(find(VSDI.condition(:,1) == code) , sel_trials);
            codemA(codi)= VSDI.condition(tricond(1),4); %mA corresponding to the code label (for the next subplot)
            roiwave = mean(waves(:,nroi,tricond),3) ;
            
            hold on;
            output = devo_peak2peak(roiwave, VSDI.timebase, window, noise, method, 0, 0);
            
            idx0= dsearchn(VSDI.timebase, 0);
            
            waveW = roiwave(idx0:output.peakidx(2));
            slopemean(nroi, codi) = mean(diff(waveW));
            
            slopemax(nroi,codi) = output.slopemax;
            peak(nroi,codi) = output.peakminusbasel;
            onset_ms(nroi,codi)= output.onsetnoise_ms;
            
            plot(VSDI.timebase,roiwave, 'color', colors(codi,:), 'linewidth', 1.3); hold on %  'linewidth', 1.8
            
            %         ylim([-0.2 .3])
            ylim(wavesylim.ave) %for high stimuli
            
            xlim([-300 600])
            ylabel('%\Delta F (trials ave)');
            
            title([ VSDI.roi.labels{nroi}])
            
            
            clear output roiwave waveW slopeval
            
            legend_labels{j} =[ num2str(VSDI.condition(tricond(1),4)) 'mA'];
            j = j+1;
            
        end %codi
        

        hold off
        
        
        ploti = ploti+1;
        
        
    end %  roi
            legend(legend_labels)

    
    sgtitle([num2str(VSDI.ref), '(', descript, ')','f:', fcode , '.', num2str(noise.SDfactor),'sd' ])
    
    name2save = fullfile(savein,['plot_D_solapadas' num2str(VSDI.ref) num2str(condi(1)) '.jpg']);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
    saveas(gcf,name2save,'jpg')
    close
end
blob()