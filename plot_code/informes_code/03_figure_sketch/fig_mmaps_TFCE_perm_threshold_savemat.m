%% s06_analysis1: FRAME-WISE STATISTICAL COMPARISON (TFCE CORRECTION)
clear

%----------------------------------------------------------------
% @SET: BASIC PARAMETERS
%----------------------------------------------------------------

ref_movie= '_12filt5' ;
% ref_movie= '_11filt4' ;
% ref_movie = '_08diff_perc_f0pre';
% ref_movie = '_09filt3';

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')

% outfield = 'peakminusbasel';% @SET: output field with choosen measure
% outfield = 'slopemax'; % mean slope from idx = 0 to peak
outfield = 'wmean';
% outfield = 'wmeanslope'; %mean slope in window (minus baseline slope)

%----------------------------------------------------------------
% @SET: FUNCTION PARAMETERS
%----------------------------------------------------------------
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

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

reject_on = [1];  %@ SET
% Subsettings:
setting.manual_reject = 0; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET+
setting.force_include = 0; %@ SET


% END OF SETTINGS


% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

% end of user_settings


%
for blocki = [4 5] %4:length(fast_condition_list)
    clearvars -except blocki fast_condition_list outfield ref_movie  savein  feedf reject_on setting path

    
    nfish = fast_condition_list{blocki,1};
    
    trial_kinds = fast_condition_list{blocki,2};
    cond_def =fast_condition_list{blocki,3};
    
    
    [VSDI] = TORus('load',nfish);
    
    temp = TORus('loadmovie',nfish,ref_movie);
    movies = temp.data(:,:,1:end-1,:);

    %----------------------------------------------------------------
    % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
    %----------------------------------------------------------------

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
    
    %----------------------------------------------------------------
    % SELECT CASES
    %----------------------------------------------------------------
    sel_trials= [];
    
    for condi = trial_kinds %to make sure that only the conditions of interest are computed)
        condtrials = makeCol(find(VSDI.condition(:,1)==condi));
        sel_trials  = [sel_trials; condtrials];
    end
    sel_trials= sort(sel_trials);
    
    if reject_on  %@ SET
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    
    
    %----------------------------------------------------------------
    % CONFIGURATION OF PARAMETERS THAT ARE SPECIFIC FOR EACH MEASURE
    %----------------------------------------------------------------
    mmaps= [];
    localclim = [];
    flagslope = 0;
    
%     switch outfield
%         
%         % ------------------------------------------------------------------------
%         case 'peakminusbasel'
%             
%             
%             localmap = jet;
%             
%             localtitle = [num2str(VSDI.ref), 'peak-b'];
%             localname = [num2str(VSDI.ref), '-peak-b'];
%             
%             
%             % ------------------------------------------------------------------------
%         case 'slopemax'
%             
%             
%             localmap = jet;
%             
%             localtitle = [num2str(VSDI.ref), 'slopemax'];
%             localname = [num2str(VSDI.ref), '-slopemax'];
%             
%             
%             % ------------------------------------------------------------------------
%         case 'onsetnoise_ms'
%             
%             
%             localmap = flipud(jet);
%             
%             localtitle = [num2str(VSDI.ref), 'onsetnoise_m_s'];
%             localname = [num2str(VSDI.ref), '-onsetnoise_ms'];
%             
%             % ------------------------------------------------------------------------
%             
%         case 'noisethresh'
%             
%             
%             localmap = jet;
%             
%             
%             localtitle = [num2str(VSDI.ref), 'peak-b '];
%             localname = [num2str(VSDI.ref), '-peak-b'];
%             
%             % ------------------------------------------------------------------------
%         case 'slopemean' % in this case,  the measure has to be computed
%             
%             j = 1;
%             flagslope = 1;
%             
%             %GET MAX AND PLOT THEM WITH THE SAME LIMITS
% %             localmax = max(abs(measureframe(:)));
%             
%             
%             % 2. CONFIGURE THE OTHER PARAMETERS
%             localmap = jet;
% %             localclim = [0 localmax];
%             
%             localtitle = [num2str(VSDI.ref), 'slopemean '];
%             localname = [num2str(VSDI.ref), '-slopemean'];
%             
%     end %switch outfield
%     
    
    %--------------------------------------
    %2. APPLY FUNCTION TO EACH TRIAL
    %--------------------------------------
    
    % APPLY THE FUNCTION (differently for slopemean and other
    % measurements)
    
    % ----------------------------------------------------------------
    for triali = makeRow(sel_trials)
        
        movtrial = squeeze(movies(:,:,:,triali));
        
        for rowi = 1:size(movtrial,1)
            for coli = 1:size(movtrial,2)
                wave = squeeze(movtrial(rowi, coli, :));
                output = devo_peak2peak(wave, VSDI.timebase, feedf.window,[], feedf.method, 0, 0);
                
                if flagslope
                    idx0= dsearchn(VSDI.timebase, 0);%get 0 index
                    waveW = wave(idx0:output.peakidx(2));
                    slopemean = mean(diff(waveW));
                    
                    mmaps(rowi, coli, triali) = slopemean;
                    
                elseif ~flagslope
                    mmaps(rowi, coli, triali) =  output.(outfield);
                end %if flagslope
                
                clear output slopemean waveW wave
                
            end %coli
        end %rowi
        
        display(triali)
        
    end %triali
    % ----------------------------------------------------------------
    
    
    %2. GET MAX-MIN TO PLOT THEM WITH THE SAME LIMITS (if we want
    %to plot the maps)
    %                     maxval= max(abs(mmaps(:)));
    %
    %         c_lim= [-mmaps mmaps];
    %         BVmap = colormap_loadBV();
    
        
    %--------------------------------------
    % PERMUTATION OF 'mmaps' (measurement map)
    %--------------------------------------
    
    j = 1;
    
    for condi = makeRow(trial_kinds)
        
        condA = condi;
        condB = force0ending(condi);
        
        if condA == condB
            continue %skips computing the control condition
        end
        
        conditions_labels{j} = [num2str(condA) 'minus' num2str(condB)];
        
        % SELECT CONDITIONS TO COMPARE AND LOAD MOVIES
        
        [idxA] = find(VSDI.condition(:,1)==condA);
        [idxB] = find(VSDI.condition(:,1)==condB);
        
        % sel_trials are trials already included
        idxA = intersect(idxA, sel_trials);
        idxB = intersect(idxB, sel_trials);
        
        
        framesA = mmaps(:,:,idxA); 
        framesB =  mmaps(:,:,idxB);
        
        meanAct(:,:,j) = mean(framesA,3); 
        % cond_def = ['code',num2str(condA),'-', num2str(VSDI.list(idxA(1)).mA),'mA''minus-control' ];
%         code_def{j} = ['code',num2str(condA),'-', num2str(condB) ];
        
        idxcond = find(VSDI.condition(:,1) == condA); idxcond = idxcond(1);
        code_def{j} = ['c=',num2str(condA),'-', num2str(condB), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
        %--------------------------------------
%         newkinds(j) = condi;

        % SETTINGS FOR PLOTTING AND SAVING
        % for saving plots
        savename = fullfile(savein,[outfield, num2str(VSDI.ref),':',code_def{j},'.jpg']); %@ SET
        
        % For the permutation test:
        % nchoosek(24, 12)
        n_perm = 1000;  %@ SET
        
        % For Plotting:
        % act_clim= [-4 4]; %@SET coloraxis of the shown colors
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % COMPUTATION:  PERMUTATION + TFCE - Diffmap with p-value based threshold
        
        % nA = length(idxA); nB = length(idxB);
        
        Data{1} = permute( framesA(:,:,:), [3 1 2]);
        Data{2} = permute( framesB(:,:,:), [3 1 2]);% control/baseline condition. Trials have to be in the first dimension
        
        % APPLY FUNCTION: PERMUTATION + TFCE
        results = ept_TFCE_VSDI(Data, 'i', n_perm); % independent trials
        
        % Maps to plot
      
        TFCEmaps.(outfield).diffmap(:,:,j) = squeeze(mean(Data{1}) - mean(Data{2})); %the first dimension is trials
        TFCEmaps.(outfield).Tobs(:,:,j) = results.Obs;
        TFCEmaps.(outfield).Pmap(:,:,j) =results.P_Values;
        TFCEmaps.(outfield).conditions = conditions_labels;
        
        TFCEmaps.dataset_ref =  ref_movie;
        TFCEmaps.perm = n_perm ;
        TFCEmaps.reject_on = reject_on;
        TFCEmaps.map.ref = VSDI.ref;
        
        clear Data results frame
        
        
        j = j+1;
    end
    
    savemat = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_TFCEperm';
    name = [num2str(VSDI.ref),'_TFCE' num2str(n_perm) 'rep_' outfield '_reject' num2str(reject_on), 'data', ref_movie ,'maps'];
    save(fullfile(savemat,name),'TFCEmaps')
    
    %% PLOT STATISTICAL DIFFERENCE BETWEEN CONDITIONS (P-THRESHOLDED)
        
    % PLOT WITH A THRESHOLD OF 0.001 IN THE second ROW
    pthresh = 0.001; %to threshold out when plotting
%     custom_map = colormap_loadBV();
    custom_map = jet;
    
    ploti = 1;
        clims= [0 0.03];

    for condi= 1:length(TFCE.codedef) 
        maps = TFCEmaps.(outfield);
        
        ax(ploti) = subplot(1, 3, ploti); %can plot only 3 activity
        backgr = VSDI.backgr(:,:,VSDI.nonanidx(1)) .* VSDI.crop.mask;
        alpha = TFCEmaps.(outfield).Pmap(:,:,ploti)< pthresh; 
        %    imagesc(imtiles(:,:,ploti))
%         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
        plot_framesoverlaid2(TFCEmaps.(outfield).diffmap(:,:,ploti),backgr, alpha  ,0, ax(ploti), clims, 1, [], custom_map);

%         idxcond = find(VSDI.condition(:,1) == newkinds(condi)); idxcond = idxcond(1);
%         tit= ['c=',num2str(newkinds(condi)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
        tit = TFCE.codedef{condi};
        set(ax(ploti),'XColor', 'none','YColor','none')
        
        ax(ploti).Title.String = tit;
        ploti = ploti+1;
    end
    sgtit = [num2str(VSDI.ref) ':' num2str(pthresh) ' p-thresh ' TFCE.map_measure '(TFCE n' TFCEmaps.perm ')'];
    sgtitle(sgtit)
    
    %             ax(9) = subplot(3,3,9)
    %             imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
    %             colormap(ax(9), bone)
    %             colorbar('off')
    %             ax(9).Visible = 'off';
    %             axis image
    
        localname = ['plot_Z_Map_stat_thresh' num2str(VSDI.ref) ' (TFCEperm' num2str(n_perm)  'rep) of_ ' outfield  '_p' num2str(pthresh) '_' ref_movie 'reject' num2str(reject_on) '.jpg'];
    
    %SAVE
    saveas(gcf, fullfile(savein,localname ), 'jpg')
    close all
    %             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
    blob()
end % blocki
blob(); pause(0.1); blob()

% end %fish loop
% ----------------------------------------------

    %% PLOT >>FROM SAVED MAT FILE<< STATISTICAL DIFFERENCE BETWEEN CONDITIONS: [P-THRESHOLDED ACTIVITY MAP]
     clear

    %--------------------------------------
    % MANUALLY load TFCEmaps structure 
    %--------------------------------------

    %--------------------------------------
    % SETTINGS 
    %--------------------------------------
    field = 'wmean';
    clims= [0 0.07];
    pthresh = 0.05; %to threshold out when plotting

    
%     custom_map = colormap_loadBV();
    custom_map = jet;
    
    
    nfish = TORus('nsubject', TFCEmaps.ref) ; 
    VSDI = TORus('load', nfish) ;

    ploti = 1;

    for condi= 1:length(TFCEmaps.(field).conditions) 
        maps = TFCEmaps.(field);
        
        ax(ploti) = subplot(1, 3, ploti); %can plot only 3 activity
        backgr = VSDI.backgr(:,:,VSDI.nonanidx(1)) .* VSDI.crop.mask;
        alpha = maps.Pmap(:,:,ploti)< pthresh; 
        
        %    imagesc(imtiles(:,:,ploti))
%         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
        plot_framesoverlaid2(maps.diffmap(:,:,ploti),backgr, alpha  ,0, ax(ploti), clims, 1, [], custom_map);

%         idxcond = find(VSDI.condition(:,1) == newkinds(condi)); idxcond = idxcond(1);
%         tit= ['c=',num2str(newkinds(condi)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
        tit = maps.conditions{condi};
        set(ax(ploti),'XColor', 'none','YColor','none')
        
        ax(ploti).Title.String = tit;
        ploti = ploti+1;
    end
    sgtit = [num2str(VSDI.ref) ':' num2str(pthresh) ' p-thresh ' field '(TFCE n' TFCEmaps.perm ')'];
    sgtitle(sgtit)
    
    %             ax(9) = subplot(3,3,9)
    %             imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
    %             colormap(ax(9), bone)
    %             colorbar('off')
    %             ax(9).Visible = 'off';
    %             axis image
    
        localname = ['plot_Z_Map_stat_thresh' num2str(TFCEmaps.ref) ' (TFCEperm' num2str(TFCEmaps.perm)  'rep) of_ ' field  '_p' num2str(pthresh) '_' TFCEmaps.dataset_ref 'reject' num2str(TFCEmaps.reject_on) '.jpg'];
    
    %SAVE
    savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
    saveas(gcf, fullfile(savein,localname ), 'jpg')
    close all
    %             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
    blob()
    
    %% PLOT >>FROM SAVED MAT FILE<< STATISTICAL DIFFERENCE BETWEEN CONDITIONS: [P-MAP]

    %--------------------------------------
    % MANUALLY load TFCEmaps structure 
    %--------------------------------------

    %--------------------------------------
    % SETTINGS 
    %--------------------------------------
    field = 'wmean';
    pthresh = 0.05; %to (maximum value to plot)

    
%     custom_map = colormap_loadBV();
    custom_map = flipud(parula);
    
    
    nfish = TORus('nsubject', TFCEmaps.ref) ; 
    VSDI = TORus('load', nfish) ;

    ploti = 1;

    for condi= 1:length(TFCEmaps.(field).conditions) 
        maps = TFCEmaps.(field);
        
        ax(ploti) = subplot(1, 3, ploti); %can plot only 3 activity
        backgr = VSDI.backgr(:,:,VSDI.nonanidx(1)) .* VSDI.crop.mask;
        alpha = maps.Pmap(:,:,ploti)< pthresh; 
        clims = [0 pthresh];
        %    imagesc(imtiles(:,:,ploti))
%         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
        plot_framesoverlaid2(maps.Pmap(:,:,ploti),backgr, alpha  ,0, ax(ploti), clims, 1, [], custom_map);

%         idxcond = find(VSDI.condition(:,1) == newkinds(condi)); idxcond = idxcond(1);
%         tit= ['c=',num2str(newkinds(condi)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
        tit = maps.conditions{condi};
        set(ax(ploti),'XColor', 'none','YColor','none')
        
        ax(ploti).Title.String = tit;
        ploti = ploti+1;
    end
    sgtit = [num2str(VSDI.ref) ':' num2str(pthresh) ' p-map ' field '(TFCE n' TFCEmaps.perm ')'];
    sgtitle(sgtit)
    
    %             ax(9) = subplot(3,3,9)
    %             imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
    %             colormap(ax(9), bone)
    %             colorbar('off')
    %             ax(9).Visible = 'off';
    %             axis image
    
        localname = ['plot_Z_Pmap' num2str(TFCEmaps.ref) ' (TFCEperm' num2str(TFCEmaps.perm)  'rep) of_ ' field  '_p' num2str(pthresh) '_' TFCEmaps.dataset_ref 'reject' num2str(TFCEmaps.reject_on) '.jpg'];
    
    %SAVE
    savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
    saveas(gcf, fullfile(savein,localname ), 'jpg')
    close all
    %             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
    blob()
    
    
%% Created: 30/08/21
% Updated: 09/09/21 
% FROM: Gent2 code 'fig_mmaps_TFCE_perm'. This is an adaptation to plot the
% measure p-thresholded