%% TEST BOXPLOT + STATS

% SETTINGS
clear
% close all

nfish = 2
    reject_on = 0;
    
    % user_settings
    path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/VSDI_ourToolbox/';
    path.TORus = '/home/tamara/Documents/MATLAB/VSDI/TORus';
    path.data = fullfile(path.TORus, 'data');
    path.list =fullfile(path.data,'BVlists');
    path.grouplist = path.TORus;
    addpath(genpath(path.rootpath));
    addpath(path.TORus);
    % end of user_settings
    
    
    %load structure + wave
    VSDI = TORus('load',nfish);
    VSDroiTS= TORus('loadwave',nfish);
    
    roikind = 'circle';
%     roikind = 'anat';
    
    cond_codes =[200:203];
    
    % PARAMETERS FOR  'peak2peak' FUNCTION
    method = 'movsum';

    % For trial-wise analysis
    window.min = [-100 100];
    window.max = [0 600];
    window.movsum = 50;%ms
    window.baseline = [-300 0];
    window.slope = 50; %ms
    
    
    noise.fr_abovenoise = 30;
    noise.SDfactor = 2;
    
    % For average-based analysis
    
    window_ave.min = [-100 100];
    window_ave.max = [0 600];
    window_ave.movsum = 50;%ms
    window_ave.baseline = [-300 0];
    window_ave.slope = 50;
    
    
    noise_ave.fr_abovenoise = 30;
    noise_ave.SDfactor = 4;%
    
    lat_limit = 1000; %for trial rejection
    
    
    %% other settings
    
    fcode = 'f309'; %@ SET
    temp = fcode(end-2:end);
    
    switch roikind
        case 'circle'
            waves = VSDroiTS.(strcat('circ_filt',temp)).data; % automatically choses the dataset according to the number of filter
        case 'anat'
            waves = VSDroiTS.(strcat('filt',temp)).data; %@ SET ...
    end
    
    
    setting.manual_reject = 0; %@ SET
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET
    setting.force_include = 1; %@ SET
    
    
    pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/boxplot+multcomp/preS_F0_filt09/multcompare_intense_noclean';
    
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
            disp(['trial:',num2str(triali)])
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
            
            measures.slopemax(nroi,triali) = output.slopemax;

            
            % get meanslope
            idx0= dsearchn(VSDI.timebase, 0);%get 0 index
            waveW = wave(idx0:output.peakidx(2));
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
    
    
    for nroi = 1:length(VSDI.roi.labels)
        roiwave = squeeze(avewave(:,nroi));
        output = devo_peak2peak(roiwave, VSDI.timebase, window_ave, noise_ave, method, 0);
        
        idx0= dsearchn(VSDI.timebase, 0);%get 0 index
        waveW = roiwave(idx0:output.peakidx(2)); %calculate from 0 to peak
        slopeval = mean(diff(waveW));
        
        avemeasures.peakminusbasel(nroi,j) = output.peakminusbasel;
        avemeasures.onsetnoise_ms(nroi,j) = output.onsetnoise_ms;
        avemeasures.meanslope(nroi,j) = slopeval;
        avemeasures.slopemax(nroi,j) = output.slopemax;
        
        clear output roiwave waveW slopeval
    end %roi
    
    % get mA corresponding to each condition (for later plotting)
    mAcond(j)=VSDI.condition(tricond(1),4);
    j = j+1;
end

%% MEASURE SELECTION

% close all
roiname1 = 'dm4m_R';
roiname2 = 'dm2_R';

% measure = {'peakminusbasel'};
% measure = {'slopemean'};
% measure = {'slopemax'};
% measure ={'onsetnoise_ms'};
measure = {'slopemax' 'slopemean' 'onsetnoise_ms' 'peakminusbasel'};

saveon = 1;


        for resulti = 1:length(measure)
            
            result= measure{resulti};
%% BOXPLOT + STAT
nroi1 = name2idx(roiname1, VSDI.roi.labels);
nroi2 = name2idx(roiname2, VSDI.roi.labels);

... reinitialize variables
    measure1 = [];
avemeasure1 = [];

measure2 = [];
avemeasure2 = [];

extratitle = '';

ccmap= lines(2);  %custom cmap

% CONFIGURATION OF PARAMETERS THAT ARE SPECIFIC FOR EACH MEASURE
switch result
    
    % ------------------------------------------------------------------------
    case 'peakminusbasel'
        
        measure1 = makeCol(squeeze(measures.peakminusbasel(nroi1,sel_trials)));
        avemeasure1 = squeeze(avemeasures.peakminusbasel(nroi1,:));
        
        measure2 = makeCol(squeeze(measures.peakminusbasel(nroi2,sel_trials)));
        avemeasure2 = squeeze(avemeasures.peakminusbasel(nroi2,:));
        
        % BUILD CONDITIONS (mA) MATRIX FOR ANOVA
        mA= [];
        mA = VSDI.condition(sel_trials,4); %mA
        
        local_mAcond = mAcond;
        
        % ------------------------------------------------------------------------
    case 'slopemean'
        
        measure1 = makeCol(squeeze(measures.meanslope(nroi1,sel_trials)));
        avemeasure1 = squeeze(avemeasures.meanslope(nroi1,:));
        
        measure2 = makeCol(squeeze(measures.meanslope(nroi2,sel_trials)));
        avemeasure2 = squeeze(avemeasures.meanslope(nroi2,:));
        
        
        % BUILD CONDITIONS (mA) MATRIX FOR ANOVA
        mA= [];
        mA = VSDI.condition(sel_trials,4); %mA
        
        local_mAcond= mAcond;

        % ------------------------------------------------------------------------
    case 'slopemax'
        
        measure1 = makeCol(squeeze(measures.slopemax(nroi1,sel_trials)));
        avemeasure1 = squeeze(avemeasures.slopemax(nroi1,:));
        
        measure2 = makeCol(squeeze(measures.slopemax(nroi2,sel_trials)));
        avemeasure2 = squeeze(avemeasures.slopemax(nroi2,:));
        
        
        % BUILD CONDITIONS (mA) MATRIX FOR ANOVA
        mA= [];
        mA = VSDI.condition(sel_trials,4); %mA
        
        local_mAcond= mAcond;

        % ------------------------------------------------------------------------
    case 'onsetnoise_ms' % leave the condition 0 out
        seltrials_onset = [];
        controltrials = find(VSDI.condition(:,1) == cond_codes(1)); %find control trials
        seltrials_onset = setdiff(sel_trials, controltrials); %leave out of the measure
        
        % BUILD CONDITIONS (mA) MATRIX FOR ANOVA
        
        mA= [];
        mA = VSDI.condition(seltrials_onset,4); %mA
        
        measure1 = makeCol(squeeze(measures.onsetnoise_ms(nroi1,seltrials_onset)));
        avemeasure1 = squeeze(avemeasures.onsetnoise_ms(nroi1,2:end)); %we skip the first condition
        
        measure2 = makeCol(squeeze(measures.onsetnoise_ms(nroi2,seltrials_onset)));
        avemeasure2 = squeeze(avemeasures.onsetnoise_ms(nroi2,2:end)); %we skip the first condition
        
        local_mAcond = mAcond(2:end);
        
end

%SCATTER PLOT (both regions together)

subplot(1,3,1)
p1= plot (avemeasure1, local_mAcond, 'color', ccmap(1,:)); hold on
p2= plot (avemeasure2, local_mAcond, '--', 'color', ccmap(2,:)); hold on
s1= scatter(avemeasure1, local_mAcond, [], ccmap(1,:), 'filled', 'displayname', roiname1); hold on
s2= scatter(avemeasure2, local_mAcond, [], ccmap(2,:), 'filled', 'displayname', roiname2);

ylabel('mA')
ylim([local_mAcond(1)*0.9 local_mAcond(end)*1.1])
xlim([min([avemeasure1 avemeasure2])*0.9 max([avemeasure1 avemeasure2])*1.1]) % xlimit according to the max min values
yticks(mAcond)
set(gca, 'ydir', 'reverse')
title('average wave')
legend([s1 s2])

%MULTICOMPARE PLOTS (one for each region)

[~, ~, stats1] = anova1(measure1,mA, 'off') ;

h = subplot(1,3,2);
[c, m , h, gnames] = multcompare_modif(stats1, 'CType', 'bonferroni')
title([VSDI.roi.labels{nroi1} 'post-hoc'])

[~, ~, stats2] = anova1(measure2,mA, 'off') ;

h2 = subplot(1,3,3);
[c, m , h2, gnames] = multcompare_modif(stats2, 'CType', 'bonferroni')
title([VSDI.roi.labels{nroi2} 'post-hoc'])


sgtitle([num2str(VSDI.ref), ':', result, '(', fcode, ')'])

% Save
if saveon
    if reject_on
        name = strcat(['2multcompare_',num2str(VSDI.ref),':',result, 'from', num2str(nroi1),'and',num2str(nroi2),'_filt', fcode,'block',num2str(cond_codes(1)), '(clean)']);
    else
        name = strcat(['2multcompare_',num2str(VSDI.ref),':',result, 'from', num2str(nroi1),'and',num2str(nroi2),'_filt', fcode,'block',num2str(cond_codes(1))]);
    end
    
    saveas(gcf,fullfile(pathsave,name),'jpg')
    display(['FIGURE SAVED: ' name])
    close
end

end %result loop

% end %result loop
%% LAST UPDATE: 17/05/21
% CREATED: 17/05/21