%% ABSOLUTE THRESHOLD ON GS rejection criterion
% If the GS is above the threshold
clear
user_settings
pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
% ------------------------------------------------------------------------
... COMPUTE GLOBAL SIGNAL FROM ALL NONAN TRIALS AND STORE THE MAX-VALUE
    % ------------------------------------------------------------------------
load(fullfile(pathmat, 'TORus_GSwaves'))

for nfish = [21]
    
    ref_movie = '_18filt6';  % '21filt6' '_18filt6'
    VSDI = TORus('load', nfish);
    VSDmov = TORus('loadmovie',nfish,ref_movie);
    movies=VSDmov.data;
    
    cropmask = VSDI.crop.mask;
    % ·······························································
    % DEPRECATED (31/03/22) from 'reject3' - %F is now calculated from the %F movies

    %     % Preallocate GS
%     GSdF = NaN(length(VSDI.timebase), VSDI.nonanidx(end));
%     maxval= NaN(length(VSDI.nonanidx(end)));
%     % Preallocate GSperc
%     GSpF = NaN(length(VSDI.timebase), length(VSDI.nonanidx));
%     maxvalperc= NaN(length(VSDI.nonanidx(end)));
%     
%     for triali = makeRow(VSDI.nonanidx)
%         % Get dF GS
%         movietrial = movies(:,:,:,triali);
%         GSdF(:,triali) = roi_TSave(movietrial,cropmask);
%         maxval(triali) = max(movmean(GSdF(:,triali), 10));
%         
%         % Get %F GS
%         F0 = movietrial(:,:,end);
%         GSpF(:,triali) = roi_TSave_percF(movietrial,cropmask, F0);
%         maxvalperc(triali) = max(movmean(GSpF(:,triali), 10));
%         
%         clear movietrial
%     end
    % ·······························································

    % Preallocate GSwave
    GSwave= NaN(length(VSDI.timebase), VSDI.nonanidx(end));
    maxval= NaN(length(VSDI.nonanidx(end)));
    
    for triali = makeRow(VSDI.nonanidx)
        % Get dF GS
        movietrial = movies(:,:,:,triali);
        twave = roi_TSave(movietrial,cropmask);
        GSwave(:,triali) = twave;
        maxval(triali) = max(movmean(twave, 10));
        
        clear movietrial twave
    end
    
    % % PLOT ALL GS WAVES WITH 'plot_waveGS'
    % name = ['GS waveplot (dF)' num2str(VSDI.ref) ref_movie];
    % plot_wavesGS(GSdF, VSDI.nonanidx, [], 15, name)
    % save_currentfig(pathplot, name)
    %
    %
    % name = ['GS waveplot (%F)' num2str(VSDI.ref) ref_movie];
    % plot_wavesGS(GSpF, VSDI.nonanidx, [], 15, name)
    % save_currentfig(pathplot, name)
    %
    
    blob()
    
    display([num2str(VSDI.ref) 'maxval in' VSDmov.units '=' num2str(max(maxval))]);
    
    % SAVE GS in structure
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    GS.(field).wave{nfish} = GSwave;
    GS.(field).max{nfish} = maxval;
    
    clear movies GSdF GSpF maxval maxvalperc
end %nfish

save(fullfile(pathmat, 'TORus_GSwaves'), 'GS')

% dim = size(GS);
% test = reshape(GS, [dim(1) 1 dim(2)]);
%  plot_wavesplot(test, 'crop', '', VSDI.nonanidx, 10)






%% 'reject4' : define movie of reference 
clear
user_settings

% -----------------------------
% DEFINE MOVIE OF REFERENCE
% -----------------------------

for nfish = [1:4 8:12]
    VSDI = TORus('load',nfish);
    VSDI.reject4.refmovie = '_21filt6';
    TORus('save', VSDI)
    clear VSDI
end


%% 'reject4' -------- CUMMULATIVE ERROR DEVIATION FROM MEAN IN GS rejection criterion -
% In this one, the GS criterion is applied second-wise and only in a
% total timewindow of 2s
% movies: _21filt6 %d

% GLOBAL SIGNAL METHOD
%  BASED ON THE AVERAGE SIGNAL OF THE WHOLE CROPPED BRAIN IN THE TRIALS
%  INCLUDED IN THE STUDY (VSDI.trials_in)
% Method by (Cheng, 2006) applied to the whole brain.

% Extract whole brain averaged (GS) from all trials from each conditions
% and average across trials
% substract from

clear
user_settings

std_factor = 2;
mode = 'check'; %'check' 'save'

pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))


% pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/';

for nfish =  21 %1:4 8:11
    
    clearvars -except nfish std_factor saveplot ref_movie GS
    
    VSDI = TORus('load',nfish);
    
    ref_movie = VSDI.reject4.refmovie; 
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    
    GSwaves = GS.(field).wave{nfish}; %takes the %F
    
    % WATCH-OUT:delete previous values
    fieldn = ['GSdeviat', num2str(std_factor),'sd']; fieldn= strrep(fieldn,'.','');
    
    VSDI.reject3.ref = ref_movie;
        
    
    %     try
    %     VSDI.reject3.(fieldn) = [];
    %     end
    %copy GS
    
    % ---------------------------------------------------------
    ... FOR BASELINE
        % ---------------------------------------------------------
    onsetidx = dsearchn(VSDI.timebase, 0);
    idx2000 = dsearchn(VSDI.timebase, 2000);
    %     idx3000 = dsearchn(VSDI.timebase, 3000);
    
    % find aberrant trials in baseline period for all trials
    
    
    all_idx = VSDI.nonanidx; % it will include all trials from all conditions
    
    [out_basel] = find_aberrant(GSwaves(1:onsetidx,:), all_idx,  std_factor, 0); title([num2str(VSDI.ref) '.Basel GS from rejected(grey)from all included trials'])
    
    % -------------------------------------------------------
    ... FOR POST-S
        % -------------------------------------------------------
    
    % find aberrant trials for each condition in the post-Stimulus period
    cond_codes = unique(VSDI.condition(:,1)); cond_codes = cond_codes(~isnan(cond_codes));
    
    out_specific = [];
    
    for ci = 1:numel(cond_codes)
        codi  = cond_codes(ci);
        idx = find(VSDI.condition(:,1) == codi);
        temp_aberrant = find_aberrant(GSwaves(onsetidx:idx2000,:), idx, std_factor, 1); title([num2str(VSDI.ref) '.Post-S GS  from rejected(grey)from cond' num2str(cond_codes(ci))])
        
        out_specific = [out_specific temp_aberrant];
    end
    
    reject_idx = sort(unique([out_basel out_specific]));
    
    criteria = ['GSdeviat', num2str(std_factor),'sd']; criteria= strrep(criteria,'.','');
    name = ['REJ criteria 4 - ' num2str(VSDI.ref) ref_movie criteria];
    
    % store in structure
    fieldn = ['GSdeviat', num2str(std_factor),'sd']; fieldn= strrep(fieldn,'.','');
    VSDI.reject4.(fieldn) = makeCol( sort(union(out_basel, out_specific)));
    VSDI.reject4.refmovie =  ref_movie; 
    
    TORus('save', VSDI)
    
    clear   aberrant_idx reject rejectB rejectP...
        aberrant_idxall aberrant_idx...
        GSwave VSDI
end



%% ----------------------------------------------------------
... 'reject4' --------- ABSOLUTE-THRESHOLD REJECT CRITERIA
% ----------------------------------------------------------
clear
user_settings
pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))

for nfish = [12]
    
    VSDI = TORus('load',nfish);
    ref_movie = VSDI.reject4.refmovie; %to avoid mistakes the first criterion computed sets the reference movie
    
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    
    GSwaves = GS.(field).wave{nfish};
    maxval = GS.(field).max{nfish};
    ref_movie = VSDI.reject4.refmovie; %to avoid mistakes the first criterion computed sets the reference movie

    
    %     abs_thr = max(max_pF)*0.75;
    abs_thr = 0.25;
    reject_idx = [];
    
    for triali = makeRow(VSDI.nonanidx)
        if maxval(triali) >= abs_thr
            reject_idx = [reject_idx; triali];
        end
%         tit= ['test' num2str(VSDI.ref) 'trial' num2str(triali)];
%         plot(GSwaves(:,triali))
%         title(tit)
%         saveas(gcf, fullfile(pathmat, [tit, '.jpg']), 'jpg')
%         close
    end
    
    reject_idx = sort(reject_idx);
    
    %     % PLOT ALL GS WAVES WITH 'plot_waveGS'
    %     criteria = 'Abs-thresh 0.75 GSmax';
    %     name = ['REJ criteria 1 - ' num2str(VSDI.ref) ref_movie criteria];
    %     plot_wavesGS(GSwaves, VSDI.nonanidx, reject_idx, 15, name)
    %     save_currentfig(pathplot, name)
    %
    %     rejectionAbs075{nfish} = reject_idx;
    VSDI.reject4.GSabs025 = makeCol(reject_idx);
    %
            TORus('save', VSDI)
    
end


% fieldname = ['GSabs',num2str(abs_thr)]; fieldname = strrep(fieldname,'.','');
% VSDI.reject2.ref = ref_movie;
% VSDI.reject2.(fieldname)= reject_idx;
% TORus('save', VSDI)
%% 'reject3' -------- CUMMULATIVE ERROR DEVIATION FROM MEAN IN GS rejection criterion -
% In this one, the GS criterion is applied second-wise and only in a
% total timewindow of 2s

% GLOBAL SIGNAL METHOD
%  BASED ON THE AVERAGE SIGNAL OF THE WHOLE CROPPED BRAIN IN THE TRIALS
%  INCLUDED IN THE STUDY (VSDI.trials_in)
% Method by (Cheng, 2006) applied to the whole brain.

% Extract whole brain averaged (GS) from all trials from each conditions
% and average across trials
% substract from

clear
user_settings

std_factor = 2;
saveplot = 1;

pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))

ref_movie = '_18filt6';


pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/';

for nfish = [21]
    
    VSDI = TORus('load',nfish);
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    
%     GSwaves = GS.(field).pF{nfish}; %takes the %F !!! BEFORE 02/09/22
    GSwaves = GS.(field).wave{nfish}; %takes the %F !!! FROM 02/09/22 ON -TAKE ORIGINAL UNITS
    
    % WATCH-OUT:delete previous values
    fieldn = ['GSdeviat', num2str(std_factor),'sd']; fieldn= strrep(fieldn,'.','');
    
    VSDI.reject3.ref = ref_movie;
        
    
    %     try
    %     VSDI.reject3.(fieldn) = [];
    %     end
    %copy GS
    
    % ---------------------------------------------------------
    ... FOR BASELINE
        % ---------------------------------------------------------
    onsetidx = dsearchn(VSDI.timebase, 0);
    idx2000 = dsearchn(VSDI.timebase, 2000);
    %     idx3000 = dsearchn(VSDI.timebase, 3000);
    
    % find aberrant trials in baseline period for all trials
    
    
    all_idx = VSDI.nonanidx; % it will include all trials from all conditions
    
    [out_basel] = find_aberrant(GSwaves(1:onsetidx,:), all_idx,  std_factor, 0); title([num2str(VSDI.ref) '.Basel GS from rejected(grey)from all included trials'])
    
    % -------------------------------------------------------
    ... FOR POST-S
    % -------------------------------------------------------
    
    % find aberrant trials for each condition in the post-Stimulus period
    cond_codes = unique(VSDI.condition(:,1)); cond_codes = cond_codes(~isnan(cond_codes));
    
    out_specific = [];
    
    for ci = 1:numel(cond_codes)
        codi  = cond_codes(ci);
        idx = find(VSDI.condition(:,1) == codi);
        temp_aberrant = find_aberrant(GSwaves(onsetidx:idx2000,:), idx, std_factor, 1); title([num2str(VSDI.ref) '.Post-S GS  from rejected(grey)from cond' num2str(cond_codes(ci))])
        
        out_specific = [out_specific temp_aberrant];
    end
    
        if saveplot == 1
    
            h =  findobj('type','figure');
            n = length(h);
    
            for ii = 1:n
                name = strcat(num2str(VSDI.ref),'GSpost_rejected', 'c', num2str(h(ii).Number));
                saveas(h(ii), fullfile(pathplot, [name, '.jpg']),'jpg')
                close(h(ii))
            end
            clear h n
    
        end
    
    reject_idx = sort(unique([out_basel out_specific]));
    
    criteria = ['GSdeviat', num2str(std_factor),'sd']; criteria= strrep(criteria,'.','');
    name = ['REJ criteria 3 - ' num2str(VSDI.ref) ref_movie criteria];
    
    % store in structure
    fieldn = ['GSdeviat', num2str(std_factor),'sd']; fieldn= strrep(fieldn,'.','');
    VSDI.reject3.(fieldn) = makeCol( sort(union(out_basel, out_specific)));
    
    TORus('save', VSDI)
    
    clear   aberrant_idx reject rejectB rejectP...
        aberrant_idxall aberrant_idx...
        GSdF VSDI
end

clearvars -except GS rejectionGSdeviation rejectionAbs075


%% ----------------------------------------------------------
... 'reject3' --------- ABSOLUTE-THRESHOLD REJECT CRITERIA
% ----------------------------------------------------------
clear
user_settings
pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))
ref_movie = '_18filt6';

for nfish =  4
    
    VSDI = TORus('load',nfish);
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    GSwaves = GS.(field).pF{nfish};
    maxval = GS.(field).maxpF{nfish};
    
    %     abs_thr = max(max_pF)*0.75;
    abs_thr = 0.0025;
    reject_idx = [];
    
    for triali = makeRow(VSDI.nonanidx)
        if maxval(triali) >= abs_thr
            reject_idx = [reject_idx; triali];
            
        end
%         tit= ['test' num2str(VSDI.ref) 'trial' num2str(triali)];
%         plot(GSwaves(:,triali))
%         title(tit)
%         saveas(gcf, fullfile(pathmat, [tit, '.jpg']), 'jpg')
%         close
    end
    
    reject_idx = sort(reject_idx);
    
    %     % PLOT ALL GS WAVES WITH 'plot_waveGS'
    %     criteria = 'Abs-thresh 0.75 GSmax';
    %     name = ['REJ criteria 1 - ' num2str(VSDI.ref) ref_movie criteria];
    %     plot_wavesGS(GSwaves, VSDI.nonanidx, reject_idx, 15, name)
    %     save_currentfig(pathplot, name)
    %
    %     rejectionAbs075{nfish} = reject_idx;
    VSDI.reject3.GSabs025 = makeCol(reject_idx);
    %
            TORus('save', VSDI)
    
end


% fieldname = ['GSabs',num2str(abs_thr)]; fieldname = strrep(fieldname,'.','');
% VSDI.reject2.ref = ref_movie;
% VSDI.reject2.(fieldname)= reject_idx;
% TORus('save', VSDI)

%% PLOT ALL REJECTED from certain conditions

% crit_1 'absolute threshold'
% crit_2 = ''

% %% PRINT REJECTED IN EXCEL
%
% clear
%
% for nfish =  [2 3 4]
%     VSDI = TORus('load',nfish);
%
%
%     set.manual_reject = 1;
%     set.GSabsthres_reject =1;
%     set.GSmethod_reject=1;
%     set.force_include =1;
%
%
%     out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject';
%     out.name = 'rejected2';
%     export.excel = 1;
%
%
%     rejectidx = [];
%
%     rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
%
%     rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];
%
%     rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];
%
%     rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
%
%
%     rejectidx = sort(unique(rejectidx));
%
%     localoutput(:,3) = rejectidx;
%     localoutput(1:end,2) =NaN;
%
%     localoutput(:,1) = VSDI.trialref(rejectidx);
%     localoutput(:,4) = NaN;
%     localoutput(:,5) = VSDI.condition(rejectidx,1); %condition
%
%     % WRITE EXCEL
%     if export.excel == 1
%         out.sheet = [num2str(VSDI.ref)];
%         excelname = fullfile(out.folder,['rejected_trials',out.name,'.xlsx']);
%         writematrix (localoutput,excelname,'sheet',out.sheet)
%     end
%     clear localoutput rejected
% end %nfish
%
%

%% LOAD BOTH REJECTION CRITERIA INDEXES AND PLOT TO COMPARE
pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))

pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/';


for nfish = [2 8:12]
    
    ref_movie = '_18filt6';
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    GSwaves = GS.(field).pF{nfish};
    
    VSDI = TORus('load',nfish);
    
    all1 = sort(unique([VSDI.reject.GSabs025; VSDI.reject.GSdeviat2sd]));
    
    name = ['reject1-' num2str(VSDI.ref)];
    plot_wavesGS(GSwaves, VSDI.nonanidx, all1, 20, name)
    save_currentfig(pathplot, name)
    
    
    all2 = sort(unique([VSDI.reject2.GSabs025; VSDI.reject3.GSdeviat18sd]));
    
    name =  ['reject2 y 3-' num2str(VSDI.ref)];
    plot_wavesGS(GSwaves, VSDI.nonanidx, all2, 20, name)
    save_currentfig(pathplot, name)
    
    clear GSwaves all1 all2
    
end

%% CONDITION-WISE PLOT OF REJECTED
pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))

pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/';

for nfish = [2 8:12]
    
    ref_movie = '_18filt6';
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    GSwaves = GS.(field).pF{nfish};
    
    VSDI = TORus('load',nfish);
    
    all3 = sort(unique([VSDI.reject3.GSabs025; VSDI.reject3.GSdeviat2sd]));
    
    name = ['reject3-' num2str(VSDI.ref)];
    plot_wavesGS(GSwaves, VSDI.nonanidx, all3, 20, name)
    save_currentfig(pathplot, name)
    
    
    clear GSwaves all1 
    
end

%% PLOT Dm4 and rejected

pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/';

for nfish = [2 8:12]
    ref_movie = '_18filt6';
    VSDI = TORus('load', nfish);
    VSDmov = TORus('loadmovie',nfish,ref_movie);
    movies=VSDmov.data;
    
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    
    % get conditions
    conditions = unique(VSDI.condition(:,1));
    conditions = conditions(~isnan(conditions));
    
    % get rejection indexes
    all1 = sort(unique([VSDI.reject.GSabs025; VSDI.reject.GSdeviat2sd]));
    all2 = sort(unique([VSDI.reject3.GSabs025; VSDI.reject3.GSdeviat18sd]));
    
    idxDm4 =name2idx('dm4m_R2', roilabels);
    cropmask = VSDI.roi.circle.mask(:,:,idxDm4);
    
    % get waves
    for condi = makeRow(conditions)
        
        seltrials = find(VSDI.condition(:,1) == condi);
        
        for triali = makeRow(seltrials)
            Y = movies(:,:,:,triali);
            F0 = VSDI.backgr(:,:,triali);
            dm4wave(:,triali) =  roi_TSave_percF_roiwise(Y,cropmask,F0);
            clear Y F0
        end
        
        % PLOT CONDITION TRIALS
        name = ['reject1-' num2str(VSDI.ref) 'cond' num2str(condi)];
        plot_wavesGS(dm4wave, seltrials, all1, 15, name)
        save_currentfig(pathplot, name)
        
        name =  ['reject 3-' num2str(VSDI.ref)  'cond' num2str(condi)];
        plot_wavesGS(dm4wave, seltrials, all2, 15, name)
        save_currentfig(pathplot, name)
        
        clear dm4wave
        
    end
    
    clear all1 all2 movies conditions
    
end

%% ----------------------------------------------------------
... reject2' -------- ABSOLUTE-THRESHOLD REJECT CRITERIA
    % ----------------------------------------------------------
pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))
ref_movie = '_18filt6';

for nfish = [1:2 8:12]
    
    VSDI = TORus('load',nfish);
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    GSwaves = GS.(field).pF{nfish};
    maxval = GS.(field).maxpF{nfish};
    
    %     abs_thr = max(max_pF)*0.75;
    abs_thr = 0.0025;
    reject_idx = [];
    
    for triali = makeRow(VSDI.nonanidx)
        if maxval(triali) >= abs_thr
            reject_idx = [reject_idx; triali];
            
        end
        tit= ['test' num2str(VSDI.ref) 'trial' num2str(triali)];
        plot(GSwaves(:,triali))
        title(tit)
        saveas(gcf, fullfile(pathmat, [tit, '.jpg']), 'jpg')
        close
    end
    
    reject_idx = sort(reject_idx);
    
    %     % PLOT ALL GS WAVES WITH 'plot_waveGS'
    %     criteria = 'Abs-thresh 0.75 GSmax';
    %     name = ['REJ criteria 1 - ' num2str(VSDI.ref) ref_movie criteria];
    %     plot_wavesGS(GSwaves, VSDI.nonanidx, reject_idx, 15, name)
    %     save_currentfig(pathplot, name)
    %
    %     rejectionAbs075{nfish} = reject_idx;
    VSDI.reject2.GSabs025 = makeCol(reject_idx);
    %
            TORus('save', VSDI)
    
end


% fieldname = ['GSabs',num2str(abs_thr)]; fieldname = strrep(fieldname,'.','');
% VSDI.reject2.ref = ref_movie;
% VSDI.reject2.(fieldname)= reject_idx;
% TORus('save', VSDI)


%% 'reject2' -------- CUMMULATIVE ERROR DEVIATION FROM MEAN IN GS rejection criterion -
% GLOBAL SIGNAL METHOD
%  BASED ON THE AVERAGE SIGNAL OF THE WHOLE CROPPED BRAIN IN THE TRIALS
%  INCLUDED IN THE STUDY (VSDI.trials_in)
% Method by (Cheng, 2006) applied to the whole brain.

% Extract whole brain averaged (GS) from all trials from each conditions
% and average across trials
% substract from

% STEP 0 Preprocess and crop brain (p10)
% STEP 1 GS for each trial. Average across trials (totalGS)

% *function to compute totalGS and total-Std.
% STEP 2 (for each trial): extract totalGS, compute the standart deviation of the residuals
% STEP 3 : compute the mean of the SD from each trial
% STEP 4 (for each trial): Discard if std of the residuals for that trial
% is > 2 times the mean of the SD
%^* in a loop

% If the GS is above the threshold

std_factor = 2;
saveplot = 1;
pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';

ref_movie = '_18filt6';

for nfish = [2 8:12]
    
    VSDI = TORus('load',nfish);
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    
    GSwaves = GS.(field).pF{nfish}; %takes the %F
    
    %
    %     % CAREFUL: DELETES PREVIOUS VALUES
    fieldn = ['GSdeviat', num2str(std_factor),'sd']; fieldn= strrep(fieldn,'.','');
    VSDI.reject2.(fieldn) = []; %delete previously calculated one
    
    
    %copy GS
    
    % ---------------------------------------------------------
    ... FOR BASELINE
        % ---------------------------------------------------------
    onsetidx = dsearchn(VSDI.timebase, 0);
    endidx = dsearchn(VSDI.timebase, 3000);
    % find aberrant trials in baseline period for all trials
    
    
    all_idx = VSDI.nonanidx; % it will include all trials from all conditions
    
    [out_basel] = find_aberrant(GSwaves(1:onsetidx,:), all_idx,  std_factor, 0); title([num2str(VSDI.ref) '.Basel GS from rejected(grey)from all included trials'])
    
    
    % -------------------------------------------------------
    ... FOR POST-S
        % -------------------------------------------------------
    
    % find aberrant trials for each condition in the post-Stimulus period
    cond_codes = unique(VSDI.condition(:,1)); cond_codes = cond_codes(~isnan(cond_codes));
    
    out_specific = [];
    
    for ci = 1:numel(cond_codes)
        codi  = cond_codes(ci);
        idx = find(VSDI.condition(:,1) == codi);
        temp_aberrant = find_aberrant(GSwaves(onsetidx:end,:), idx, std_factor, 0); title([num2str(VSDI.ref) '.Post-S GS  from rejected(grey)from cond' num2str(cond_codes(ci))])
        out_specific = [out_specific temp_aberrant];
    end
    
    %     if saveplot == 1
    %
    %         h =  findobj('type','figure');
    %         n = length(h);
    %
    %         for ii = 1:n
    %             name = strcat(num2str(VSDI.ref),'GSpost_rejected', 'c', num2str(h(ii).Number));
    %             saveas(h(ii), fullfile(pathplot, [name, '.jpg']),'jpg')
    %             close(h(ii))
    %         end
    %         clear h n
    %
    %     end
    
    reject_idx = sort(unique([out_basel out_specific]));
    % criteria = ['GSdeviat_BASELINE', num2str(std_factor),'sd']; criteria= strrep(fieldn,'.','');
    % name = ['REJ criteria 2' num2str(VSDI.ref) ref_movie criteria];
    % plot_wavesGS(GSwaves, VSDI.nonanidx, out_basel, 15, name)
    % save_currentfig(pathplot, name)
    %
    % criteria = ['GSdeviat_SPECIFIC', num2str(std_factor),'sd']; criteria= strrep(fieldn,'.','');
    % name = ['REJ criteria 2' num2str(VSDI.ref) ref_movie criteria];
    %plot_wavesGS(GSwaves, VSDI.nonanidx, out_basel, 15, name)
    % save_currentfig(pathplot, name)
    
    criteria = ['GSdeviat', num2str(std_factor),'sd']; criteria= strrep(criteria,'.','');
    name = ['REJ criteria 2 - ' num2str(VSDI.ref) ref_movie criteria];
    %     plot_wavesGS(GSwaves, VSDI.nonanidx, reject_idx, 15, name)
    %     save_currentfig(pathplot, name)
    %
    
    % store in structure
    fieldn = ['GSdeviat', num2str(std_factor),'sd']; fieldn= strrep(fieldn,'.','');
    VSDI.reject2.(fieldn) = makeCol( sort(union(out_basel, out_specific)));
    
    TORus('save', VSDI)
    
    clear   aberrant_idx reject rejectB rejectP...
        aberrant_idxall aberrant_idx...
        GSdF VSDI
end

clearvars -except GS rejectionGSdeviation rejectionAbs075
%% LOAD BOTH REJECTION CRITERIA INDEXES AND PLOT TO COMPARE
pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))

pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/';


for nfish = [2 8:12]
    
    ref_movie = '_18filt6';
    field = [ref_movie(4:end) '_' ref_movie(2:3)];
    GSwaves = GS.(field).pF{nfish};
    
    VSDI = TORus('load',nfish);
    
    all1 = sort(unique([VSDI.reject.GSabs025; VSDI.reject.GSdeviat2sd]));
    
    name = ['reject1-' num2str(VSDI.ref)];
    plot_wavesGS(GSwaves, VSDI.nonanidx, all1, 20, name)
    save_currentfig(pathplot, name)
    
    
    all2 = sort(unique([VSDI.reject2.GSabs025; VSDI.reject2.GSdeviat15sd]));
    
    name =  ['reject2-' num2str(VSDI.ref)];
    plot_wavesGS(GSwaves, VSDI.nonanidx, all2, 20, name)
    save_currentfig(pathplot, name)
    
    clear GSwaves all1 all2
    
end
%% Update history:
% 31/03/22 : reject4
% -/-/- : reject3
% Created: 28/01/2022
