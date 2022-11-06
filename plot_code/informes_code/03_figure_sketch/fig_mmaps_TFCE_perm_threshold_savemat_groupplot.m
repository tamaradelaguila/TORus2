%% s06_analysis1: FRAME-WISE STATISTICAL COMPARISON (TFCE CORRECTION)


clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)


% for blocki = 1:length(groupplot)
for nfish = [27:30] %29 sólo 701 

 clearvars -except nfish path       
%----------------------------------------------------------------
% @SET: BASIC PARAMETERS
%----------------------------------------------------------------
ref_movie= '_18filt6' ; % '_16diff_f0pre' ; '_20perc_f0pre'; '_18filt5'

% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new.mat') % ATTENTION: select cases manually (there are repetitions)

reject_on = [3];
trial_kinds = [0 701 ]; % blank conditions needs to be included and be the first for the code to properly work

% outfield = 'peakminusbasel';% @SET: output field with choosen measure
% outfield = 'slopemax'; % mean slope from idx = 0 to peak
outfield = 'wmean';
% outfield = 'wslope'; %mean slope in window (minus baseline slope)

if strcmpi(outfield, 'wslope')
flagslope = 1;
else
    flagslope = 0;
end

%----------------------------------------------------------------
% @SET: FUNCTION PARAMETERS
%----------------------------------------------------------------
feedf.window.min = [-100 100]; % 'feed-function' structure
% feedf.window.max = [0 600]; %muted 27/10/22
feedf.window.max = [0 800];

feedf.window.movsum = 50;
feedf.window.basel = [-25 0];
feedf.window.slope=50;
feedf.window.wmean=[0 350];

slope.window = [0, 200]; %ms
% feedf.noise.fr_abovenoise = 30;
% feedf.noise.SDfactor = 2;
% feedf.noise.basel = [-200 0]; %it won't be used by the function, but will be used to manually get the baseline

% Window For average-based analysis

% feedf.window_ave = feedf.window;
%
% feedf.noise_ave = feedf.noise;
% feedf.noise_ave.SDfactor = 4;% SET differences

feedf.method = 'movsum';

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

% Subsettings:
setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 0; %@ SET+
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

% savemat = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_TFCEperm/';
savemat = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/2_pmaps';
        
        %----------------------------------------------------------------
        % START CODE
        %----------------------------------------------------------------
        VSDI = TORus('load', nfish) ;

        % SETTINGS THAT ALLOW LOADING STRUCTURE
        n_perm = 1000;  %@ SET
        
        
        %----------------------------------------------------------------
        % LOAD THE STRUCTURE WITH PERMUTATION RESULTS (if exists)
        %----------------------------------------------------------------
        
        name = [num2str(VSDI.ref) '_TFCE' num2str(n_perm) 'rep_reject' num2str(reject_on), 'data', ref_movie ,'maps.mat'];
        
        try
            load(fullfile(savemat,name))
        catch
            disp(['A new result structure will be created and saved as mat-file as: ' name]);
        end
        
        %         trial_kinds = [401:404];
%         trial_kinds = unique(VSDI.condition(:,1));
%         trial_kinds= trial_kinds( ~isnan(trial_kinds));

        temp = TORus('loadmovie',nfish,ref_movie);
        movies = temp.data(:,:,1:end-1,:);
        
        %----------------------------------------------------------------
        % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
        %----------------------------------------------------------------
        reject_on=3;
        rej = 'reject' ;
        if reject_on > 1
            rej = [rej num2str(reject_on)];
        end


        rejectidx = [];
        
        if setting.manual_reject
            rejectidx = [rejectidx  makeRow(VSDI.(rej).manual)];
        end
        
        if setting.GSabsthres_reject
            rejectidx = [rejectidx  makeRow(VSDI.(rej).GSabs025)];
            
        end
        
        if setting.GSmethod_reject
            rejectidx = [rejectidx makeRow(VSDI.(rej).GSdeviat2sd)];
        end
        
        if setting.force_include
            rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
            
        end
        
        rejectidx = sort(unique(rejectidx));
        
        %----------------------------------------------------------------
        % SELECT CASES
        %----------------------------------------------------------------
        sel_trials= [];
        
        for condi = makeRow(trial_kinds) %to make sure that only the conditions of interest are computed)
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
        
        slope.windowidx = dsearchn(VSDI.timebase, makeCol(slope.window));
        slope.windowidx = [slope.windowidx(1) slope.windowidx(end)]; 

        
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
%                         idx0= dsearchn(VSDI.timebase, 0);%get 0 index
%                         idxend = output.peakidx(2);
                        idx0 = slope.windowidx(1);
                        idxend = slope.windowidx(end);
                        waveW = wave(idx0:idxend);
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
        
        condB = trial_kinds(1);
        for condi = makeRow(trial_kinds(2:end))
            
            condA = condi;
%             condB = force0ending(condi); % when the block has its corresponding blank

            
            if condA == condB
                continue %skips computing the control condition
            end
            
            conditions_labels{j} = [num2str(condA) 'minus' num2str(condB)];
            conditions_number(j) = condA;
            
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
%             savename = fullfile(savein,[outfield, num2str(VSDI.ref),'rej', num2str(reject_on), ':',code_def{j},'.jpg']); %@ SET
            
            % For the permutation test:
            % nchoosek(24, 12)
            
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
            TFCEmaps.map.ref = VSDI.ref;
            
            TFCEmaps.(outfield).diffmap(:,:,j) = squeeze(mean(Data{1}) - mean(Data{2})); %the first dimension is trials
            TFCEmaps.(outfield).Tobs(:,:,j) = results.Obs;
            TFCEmaps.(outfield).Pmap(:,:,j) =results.P_Values;
            TFCEmaps.(outfield).meanmap(:,:,j) = squeeze(mean(Data{1})); %TO DO!!!!
            TFCEmaps.(outfield).conditions = conditions_labels;
            TFCEmaps.(outfield).conditions_number = conditions_number;
            
            TFCEmaps.dataset_ref =  ref_movie;
            TFCEmaps.perm = n_perm ;
            TFCEmaps.reject_on = reject_on;
            
            clear Data results frame
            
            j = j+1;
        end %condi
        
        name = [num2str(VSDI.ref) '_TFCE' num2str(n_perm) 'rep_reject' num2str(reject_on), 'data', ref_movie ,'maps'];
        save(fullfile(savemat,name),'TFCEmaps')
        
        %% PLOT STATISTICAL DIFFERENCE BETWEEN CONDITIONS (P-THRESHOLDED)
        
        %     % PLOT WITH A THRESHOLD OF 0.001 IN THE second ROW
        %     pthresh = 0.05; %to threshold out when plotting
        % %     custom_map = colormap_loadBV();
        %     custom_map = jet;
        %
        %         clims= [];
        %
        %     for condi= 1:length(TFCEmaps.codedef)
        %         maps = TFCEmaps.(outfield);
        %
        %         ax(1) = subplot(1, 2, 1); %can plot only 3 activity
        %         backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));
        %         alpha = TFCEmaps.(outfield).Pmap(:,:,ploti)< pthresh;
        %         %    imagesc(imtiles(:,:,ploti))
        % %         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
        %         plot_framesoverlaid2(TFCEmaps.(outfield).diffmap(:,:,ploti),backgr, alpha  ,0, ax(ploti), clims, 1, [], custom_map);
        %
        % %         idxcond = find(VSDI.condition(:,1) == newkinds(condi)); idxcond = idxcond(1);
        % %         tit= ['c=',num2str(newkinds(condi)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
        %         tit = TFCE.codedef{condi};
        %         set(ax(ploti),'XColor', 'none','YColor','none')
        %
        %         ax(1).Title.String = tit;
        %         ax(1).Visible = 'off';
        %
        %         ax(2) = subplot(1,2,2);
        %         map = maps.meanmap(:,:,condi);
        %
        %         imagesc(map); axis image
        %         title([outfield ' mean map'])
        %         ax(2).Visible = 'off';
        %
        %     end
        %     sgtit = [num2str(VSDI.ref) ':' num2str(pthresh) ' p-thresh ' TFCE.map_measure '(TFCE n' TFCEmaps.perm ') rej' num2str(reject_on)];
        %     sgtitle(sgtit)
        %
        %     %             ax(9) = subplot(3,3,9)
        %     %             imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
        %     %             colormap(ax(9), bone)
        %     %             colorbar('off')
        %     %             ax(9).Visible = 'off';
        %     %             axis image
        %
        %         localname = ['P_Map_thresh' num2str(VSDI.ref) 'block' num2str(condB)  ' (TFCEperm' num2str(n_perm)  'rep) of_ ' outfield  '_p' num2str(pthresh) '_' ref_movie 'reject' num2str(reject_on) '.jpg'];
        %
        %     %SAVE
        %     saveas(gcf, fullfile(savein,localname ), 'jpg')
        %     close all
        %     %             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
        %     blob()
        pause(60*3)
    end % nfishi
% blob(); pause(0.1); blob()

% end %fish loop
% ----------------------------------------------


%% PLOT >>FROM SAVED MAT FILE<< STATISTICAL DIFFERENCE BETWEEN CONDITIONS: [P-THRESHOLDED ACTIVITY MAP]
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus/'
user_settings
cd(W)
%--------------------------------------
% MANUALLY load TFCEmaps structure
%--------------------------------------

% nfish = TORus('who', TFCEmaps.ref) ; %needed if manually loaded

%--------------------------------------
% SETTINGS TO LOAD DATA
%--------------------------------------
% outfield = 'wmean';

n_perm = 1000;
clims= [];
% outfield = 'peakminusbasel';
% outfield = 'wmean';
outfield = 'wslope';

% reject_on = 3;

%     custom_map = colormap_loadBV();
ccmap = jet;
% ccmapP = flipud(hot);
ccmapP = hot;


for nfish = [29] 
    VSDI = TORus('load', nfish) ; 
    
    for pthresh = [0.001 0.05] %to (maximum value to plot)
        
        reject_on = 3;
            for ref_movie = {'_16diff_f0pre'} %{ '_17filt5', '_18filt6'}
                
%                 matpath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_TFCEperm/';
 matpath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/2_pmaps'    ;            

                name = [num2str(VSDI.ref) '_TFCE' num2str(n_perm) 'rep_reject' num2str(reject_on), 'data', ref_movie{1} ,'maps.mat'];
                
                try
                    load(fullfile(matpath,name))
                catch
                    disp(['the file "' name '" is not found']);
                    return
                end
                
                maps = TFCEmaps.(outfield);
                
                for condi= 1:length(TFCEmaps.(outfield).conditions_number)
                    
                    figure('units','normalized','outerposition',[0 0 1 1])
                    
                    % PLOT MEANMAPS P-THRESHOLDED
                    ax(1) = subplot(1,3,1);
                    backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));% .* VSDI.crop.mask;
                    alpha = maps.Pmap(:,:,condi)< pthresh;
                    
                    clims = [];
                    %    imagesc(imtiles(:,:,ploti))
                    %         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
                    plot_framesoverlaid2(maps.diffmap(:,:,condi),backgr, alpha  ,0, ax(1), clims, 1, [], ccmap);
                    

                    ax(1).Visible = 'off';
                    
                    % PLOT P-MAP
                    ax(2) = subplot(1,3,2);
                    ax(2).Visible = 'off';
                    
                    clims = [0 pthresh];
                    %    imagesc(imtiles(:,:,ploti))
                    %         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
                    plot_framesoverlaid2(maps.Pmap(:,:,condi),backgr, alpha  ,0, ax(2), clims, 1, [], ccmapP);
                    
                    set(ax(2),'XColor', 'none','YColor','none')
                    
                    % PLOT BACKGROUND
                    ax(3) = subplot(1,3,3);
                    imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
                    axis image
                    colormap(ax(3), bone)
                    colorbar('off')
                    
                    ax(3).Visible = 'on';
                    set(ax(3),'XColor', 'none','YColor','none')
                    
                                        
                    %         idxcond = find(VSDI.condition(:,1) == newkinds(condi)); idxcond = idxcond(1);
                    %         tit= ['c=',num2str(newkinds(condi)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
                    tit = [num2str(VSDI.ref) '-' maps.conditions{condi} ': (1)' outfield ' p-thresh (2)pmap (3) backgr'];
                    
                    sgtitle(tit);
                    
                    
                    localname = ['plot_MeanMap_Pthresh_' num2str(VSDI.ref) '-' maps.conditions{condi} ' (TFCEperm' num2str(TFCEmaps.perm)  'rep) of_ ' outfield  '_p' num2str(pthresh) '_' TFCEmaps.dataset_ref 'reject' num2str(TFCEmaps.reject_on) '.jpg'];
                    %SAVE
                    savein = fullfile('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_TFCEperm/plots',outfield) ;%@ SET
                    
                    try
                        saveas(gcf, fullfile(savein,localname ), 'jpg')
                        
                    catch
                        mkdir(savein)
                        saveas(gcf, fullfile(savein,localname), 'jpg')
                    end
                    
%                     close
                    
                end % for condi
                
            end % for ref_movie
    end % for pthresh
    
end % for nfish
%             ax(9).Visible = 'off';
%             axis image

%             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
% blob()


%% PLOT >>FROM SAVED MAT FILE<< ONLY P-MAP --- HIGH QUALITY
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus/'
user_settings
cd(W)
%--------------------------------------
% MANUALLY load TFCEmaps structure
%--------------------------------------

% nfish = TORus('who', TFCEmaps.ref) ; %needed if manually loaded

%--------------------------------------
% SETTINGS TO LOAD DATA
%--------------------------------------
% outfield = 'wmean';

n_perm = 1000;
clims= [];
% outfield = 'peakminusbasel';
outfield = 'wmean';
% outfield = 'wslope';

% reject_on = 3;

%     custom_map = colormap_loadBV();
ccmap = jet;
ccmapP = flipud(hot);
% ccmapP = hot;


for nfish = [27:30] 
    VSDI = TORus('load', nfish) ; 
    
    for pthresh = [0.05 0.001] %to (maximum value to plot)
        
        reject_on = 3;
            for ref_movie = {'_16diff_f0pre'} %{ '_17filt5', '_18filt6'}  {'_20diff_f0pre'}
                
                matpath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/2_pmaps';
                
                name = [num2str(VSDI.ref) '_TFCE' num2str(n_perm) 'rep_reject' num2str(reject_on), 'data', ref_movie{1} ,'maps.mat'];
                
                try
                    load(fullfile(matpath,name))
                catch
                    disp(['the file "' name '" is not found']);
                    return
                end
                
                maps = TFCEmaps.(outfield);
                    backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));% .* VSDI.crop.mask;
                                    backgr = interp2(single(backgr),5, 'cubic');

                for condi= 1:length(TFCEmaps.(outfield).conditions_number)
                    
                    figure('units','normalized','outerposition',[0 0 1 1])
                    
                    % PLOT P-MAP
                    ax(2) = subplot(1,1,1);
                    ax(2).Visible = 'off';
                    
                    clims = [0 pthresh];
                    %    imagesc(imtiles(:,:,ploti))
                    %         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
                    alpha = maps.Pmap(:,:,condi)< pthresh;
                    alpha =  interp2(alpha,5, 'nearest');
                    act = interp2(maps.Pmap(:,:,condi),5, 'nearest');

                    plot_framesoverlaid2(act, backgr, alpha  ,0, ax(2), clims, 1, [], ccmapP);
                    
                    set(ax(2),'XColor', 'none','YColor','none')
                    
                    tit = [num2str(VSDI.ref) '-' maps.conditions{condi} ': (1)' outfield ' p-thresh (2)pmap (3) backgr'];
                    sgtitle(tit);
                    
                    localname = ['HIdef_Pmap_' num2str(VSDI.ref) '-' maps.conditions{condi} ' (TFCEperm' num2str(TFCEmaps.perm)  'rep) of_ ' outfield  '_p' num2str(pthresh) '_' TFCEmaps.dataset_ref 'reject' num2str(TFCEmaps.reject_on) '.jpg'];
                    %SAVE
%                     savein = fullfile('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_TFCEperm/plots',outfield) ;%@ SET
                    savein= fullfile(matpath, outfield);
                    
                    try
                        saveas(gcf, fullfile(savein,localname ), 'jpg')
                        
                    catch
                        mkdir(savein)
                        saveas(gcf, fullfile(savein,localname), 'jpg')
                    end
                    
                    close
                    
                end % for condi
                
            end % for ref_movie
    end % for pthresh
    
end % for nfish
%             ax(9).Visible = 'off';
%             axis image

%             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
% blob()
% %% PLOT >>FROM SAVED MAT FILE<< STATISTICAL DIFFERENCE BETWEEN CONDITIONS: [P-MAP]
%
% %--------------------------------------
% % MANUALLY load TFCEmaps structure
% %--------------------------------------
%
% nfish = TORus('who', TFCEmaps.ref) ;
% VSDI = TORus('load', nfish) ;
%
% %     custom_map = colormap_loadBV();
% custom_map = flipud(parula);
%
%
% %--------------------------------------
% % SETTINGS TO LOAD DATA
% %--------------------------------------
% field = 'wmean';
% pthresh = 0.05; %to (maximum value to plot)
% n_perm = 1000;
%
% reject_on = 0;
% ref_movie = '_17filt5';
%
% matpath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_TFCEperm';
% name = [num2str(VSDI.ref)  '_TFCE' num2str(n_perm) 'rep_' outfield '_reject' num2str(reject_on), 'data', ref_movie ,'maps.mat'];
% load(fullfile(matpath,name))
%
%
% for condi= 1:length(TFCEmaps.(field).conditions)
%     maps = TFCEmaps.(field);
%
%     backgr = VSDI.backgr(:,:,VSDI.nonanidx(1)) .* VSDI.crop.mask;
%     alpha = maps.Pmap(:,:,ploti)< pthresh;
%     clims = [0 pthresh];
%     %    imagesc(imtiles(:,:,ploti))
%     %         plot_logpmap_overlaid(maps.Pmap(:,:,condi),VSDI.backgr(:,:,VSDI.nonanidx(1)),pthresh ,0, ax(ploti), 1, flipud(parula));
%     plot_framesoverlaid2(maps.Pmap(:,:,ploti),backgr, alpha  ,0, ax(ploti), clims, 1, [], custom_map);
%
%     %         idxcond = find(VSDI.condition(:,1) == newkinds(condi)); idxcond = idxcond(1);
%     %         tit= ['c=',num2str(newkinds(condi)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
%     tit = maps.conditions{condi};
%     set(ax(ploti),'XColor', 'none','YColor','none')
%
%     ax(ploti).Title.String = tit;
%     localname = ['plot_Pmap' num2str(TFCEmaps.ref) tit ' (TFCEperm' num2str(TFCEmaps.perm)  'rep) of_ ' field  '_p' num2str(pthresh) '_' TFCEmaps.dataset_ref 'reject' num2str(TFCEmaps.reject_on) '.jpg'];
%
%     %SAVE
%
%     savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
%     saveas(gcf, fullfile(savein,localname ), 'jpg')
%     close all
%
% end %condi
% %     sgtit = [num2str(VSDI.ref) ':' num2str(pthresh) ' p-map ' field '(TFCE n' TFCEmaps.perm ')'];
% %     sgtitle(sgtit)
%
% %             ax(9) = subplot(3,3,9)
% %             imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
% %             colormap(ax(9), bone)
% %             colorbar('off')
% %             ax(9).Visible = 'off';
% %             axis image
%
% %             clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
% blob()


%% Created: 30/08/21
% Updated: 22/12/21
% Updated: 28/12/21 - to save each p-thresholded map individually
% FROM: Gent2 code 'fig_mmaps_TFCE_perm'. This is an adaptation to plot the
% measure p-thresholded