
clear
ref_wave = 'circ_filt309';
ref_movie= '_09filt3' ;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_anova' ;%@ SET
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')

outfield = 'wmean';
% 'peakminusbasel' , 'meanslope', 'wmean'

% COMMON SETTINGS FOR THE FUNCTION
feedf.window.min = [-100 100]; % 'feed-function' structure
feedf.window.max = [0 1000];
feedf.window.movsum = 50;
feedf.window.basel = [-100 0];
feedf.window.slope=50;
feedf.window.wmean=[0 250];

feedf.noise.fr_abovenoise = 30;
feedf.noise.SDfactor = 2;

% Window For average-based analysis

feedf.window_ave = feedf.window;

feedf.noise_ave = feedf.noise;
feedf.noise_ave.SDfactor = 4;% SET differences


feedf.method = 'movsum';

% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));
% END OF USER SETTINGS

% Input control
if strcmpi(outfield, 'slopemean')
    flagslope = 1;
else
    flagslope = 0;
end

%%

selroinames = {'dm4_R','dm4m_R', 'dm2_R', 'dm3_R', 'dldm_R', 'dm1_R'}; % @SET

for  block = 4%1:length(fast_condition_list)
    
    nfish = fast_condition_list{block,1};
    
    trial_kinds = fast_condition_list{block,2};
    cond_def =fast_condition_list{block,3};
    
    VSDI = TORus('load',nfish);
    %     spike = TORus('loadspike', nfish); % ECG
    VSDroiTS =TORus('loadwave',nfish);
    
    selroi =name2idx(selroinames, VSDI.roi.labels_circ); %idx of selected subset
    waves = VSDroiTS.(ref_wave).data(:,selroi,:); %@ SET
    
    
    %     output = {'trial', 'spiketime', 'cond','mA','', '%brady(count)','', 'IBIbasel', 'ibi0', 'ibi1','ibi2','ibi3','', 'estimation','dm4','dm4m','dm2','dldm'}; %header first
    
    % ----------------------------------------------
    % BRAINVISION/WAVE-BASED MEASURES
    % ----------------------------------------------
    for triali =  makeRow (VSDI.nonanidx)
        
        for nroi = 1:size(waves,2)
            wave = squeeze(waves(:, nroi,triali));
            local_output = devo_peak2peak(wave, VSDI.timebase, feedf.window,[], feedf.method, 0);
            
            
            if flagslope
                idx0= dsearchn(VSDI.timebase, 0);%get 0 index
                waveW = wave(idx0:local_output.peakidx(2));
                slopemean = mean(diff(waveW));
                msr(triali, nroi) = slopemean;
                
            elseif ~flagslope
                msr(triali, nroi) =  local_output.(outfield);
            end %if flagslope
            
            clear wave local_output
        end
        
    end %triali
    
    
    % ----------------------------------------------
    % MEAN VALUES FOR EACH CONDITION
    % ----------------------------------------------
    
    % Find the minimum number of conditions
    j = 1;
    for condi = makeRow( trial_kinds)
        sel_trials = find(VSDI.condition(:,1) == condi);
        amount(j) = length(sel_trials);
    end
    rept = min(amount); % nºof repetitions (measures of each group and measurement)
    n = rept-1;
    
    % HEADER
    coli = 1; %for conditions (levels of factor 'mA')
    for condi = makeRow(trial_kinds)
        sel_trials = find(VSDI.condition(:,1) == condi);
        sel_trials = sel_trials(1:rept) ;
        rowi = 1; % for trials (all trial from each roi subsequently)
        for roii = 1:length(selroi)
            data(rowi:rowi+n,coli) = msr(sel_trials,roii);
            rowi = rowi +rept;
            
        end % roii
        
        mA(coli) = VSDI.condition(sel_trials(1),4) ;
        
        coli = coli+1;
        clear sel_trials
    end % condi
    
    
    % ADAPT 'mA' (for x-axis labels) if there is a tone
    
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
    
    
    
    if flag_adapt
        xticklabels(labels)
    end
    
    
    [p , tb1 , stats] = anova2(data,rept);
    c = multcompare(stats)
    c = multcompare(stats,'Estimate','row')
    c = multcompare(stats,'Estimate','col')
    
    
    
    % ----------------------------------------------
    % WRITE EXCEL
    % ----------------------------------------------
    %     if export.excel == 1
    %         % write output (new sheet for each fish
    %         localoutput = cell2table(localoutput);
    %         writetable (localoutput, excelname, 'sheet', num2str(VSDI.ref))
    %     end
    
    %     name2save = fullfile(savein,['plot_C_Amplitude_trialwise' num2str(VSDI.ref) num2str(condi(1)) '.jpg']);
    %     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
    %     saveas(gcf,name2save,'jpg')
    %
    %     close
    %     clear mA meanAct stdAct  msr
    
end %nfish

blob()

%% R format

selroinames = {'dm4_R','dm4m_R', 'dm2_R', 'dm3_R', 'dldm_R', 'dm1_R'}; % @SET

for  block = 4%1:length(fast_condition_list)
    
    nfish = fast_condition_list{block,1};
    
    trial_kinds = fast_condition_list{block,2};
    cond_def =fast_condition_list{block,3};
    
    VSDI = TORus('load',nfish);
    %     spike = TORus('loadspike', nfish); % ECG
    VSDroiTS =TORus('loadwave',nfish);
    
    selroi =name2idx(selroinames, VSDI.roi.labels_circ); %idx of selected subset
    waves = VSDroiTS.(ref_wave).data(:,selroi,:); %@ SET
    
    
    %     output = {'trial', 'spiketime', 'cond','mA','', '%brady(count)','', 'IBIbasel', 'ibi0', 'ibi1','ibi2','ibi3','', 'estimation','dm4','dm4m','dm2','dldm'}; %header first
    
    % ----------------------------------------------
    % BRAINVISION/WAVE-BASED MEASURES
    % ----------------------------------------------
    for triali =  makeRow (VSDI.nonanidx)
        
        for nroi = 1:size(waves,2)
            wave = squeeze(waves(:, nroi,triali));
            local_output = devo_peak2peak(wave, VSDI.timebase, feedf.window,[], feedf.method, 0);
            
            
            if flagslope
                idx0= dsearchn(VSDI.timebase, 0);%get 0 index
                waveW = wave(idx0:local_output.peakidx(2));
                slopemean = mean(diff(waveW));
                msr(triali, nroi) = slopemean;
                
                
            elseif ~flagslope
                msr(triali, nroi) =  local_output.(outfield);
            end %if flagslope
            
            clear wave local_output
        end
        
    end %triali
    
    
    
    % ----------------------------------------------
    % MEAN VALUES FOR EACH CONDITION
    % ----------------------------------------------
    
    % Find the minimum number of conditions
    j = 1;
    intrials = [];
    for condi = makeRow( trial_kinds)
        sel_trials = find(VSDI.condition(:,1) == condi);
        amount(j) = length(sel_trials);
        
        mA(condi) = VSDI.condition(sel_trials(1),4) ;
        intrials = [intrials; sel_trials]
        clear sel_trials
    end
    intrials = sort(intrials);
    
    rept = min(amount); % nºof repetitions (measures of each group and measurement)
    n = rept-1;
    
    % HEADER
    j = 1;
    for triali = makeRow(intrials) %only include from the selected blocks
        
        for roii = 1:length(selroi)
            data_long{j,1} = msr(triali,roii);
            data_long{j,2} = VSDI.condition(triali,4);
            data_long{j,3} = roii;
            data_long{j,4} = selroinames{roii};
            
            j = j+1;
        end % roii
        
    end % condi
    
    
    % ADAPT 'mA' (for x-axis labels) if there is a tone
    
    % ----------------------------------------------
    % WRITE EXCEL
    % ----------------------------------------------
    excelname = fullfile(savein, ['long_format_forR_' outfield 'data' ref_wave '.xls']);
    csvname = fullfile(savein, ['long_format_forR_' outfield num2str(VSDI.ref) '_' num2str(trial_kinds(1))  'data' ref_wave '.csv']);
    % write output (new sheet for each fish
    data_long = cell2table(data_long);
    writetable (data_long, excelname, 'Sheet', [num2str(VSDI.ref) ])
    writetable (data_long, csvname)
    
    %     name2save = fullfile(savein,['plot_C_Amplitude_trialwise' num2str(VSDI.ref) num2str(condi(1)) '.jpg']);
    %     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
    %     saveas(gcf,name2save,'jpg')
    %
    %     close
    %     clear mA meanAct stdAct  msr
    
end %nfish