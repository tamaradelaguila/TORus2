%% CODE: TILES + WAVEs (all roi for each condition)

clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

% TO LOOP FROM STRUCTURE
...................................
    
% FOR NEW BOXPLOTS (no-pain vs pain)
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot8_electrical_intersect.mat') % ATTENTION: select cases manually (there are repetitions)

% ...FOR 'CHORRITO' BOXPLOT:
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot7_chorrito_intersect.mat') % ATTENTION: select cases manually (there are repetitions)

selroinames = {'dm4m_R2',  'dm2_R2'};%dm3 ORIGINAL

refroiname = 'dm4m_R2'; % respectb which the start/end time of the tiles will be set
% ....................................

plottiles = 1; %also plots early-peak
savetiles = 1;

plot_contour = 0;
save_contour = 0;

highdefinition = 1;

contour_of = 'peak'; % 'early' 'peak'

plot_earlypeak = 0;
save_earlypeak =0;

plotwaves=1;
savewaves =1;

keep_onlyNtrials = 0; %if >0, but keep the first 'keep_onlyNtrials' selected frames for the average

movie_ref = '_18filt6'; % input movie '_17filt5'
% movie_ref = '_17filt5'; % input movie '_17filt5'
% movie_ref = '_15filt5'; % usado para 220608

% savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/def_figs/tiles/'; %
savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/'; %

% FOR TILES AND EARLY-PEAK FRAMES
%---------------------------------------------------------------
% fact_thresh =  0.4; % @SET : limits parameters
fact_clim= 1.3;

thresh_mode = 'wavebased_thresh_max'; % 'moviebased_thresh_max' 'wavebased_thresh_max' 'wavebased_thresh_local' 'manual'. 'moviebased_thresh_max' is the one we have been using for tiles
timerange_mode = 'auto_fixed100%' ; %'auto_localwave', 'manual', 'auto_fixed', 'auto_fixed75%', 'auto_fixed100%'

analysisref = 'new4_group19_dolor_n5'; % 'new4_group19_dolor_n5' 'new4_group16_chorrito_n3'
% auto_localwave - sets the timerange according to the condition's
% reference wave

% auto_fixed - when there several conditions, it calculates the reference wave of
% all conditions, takes the highest wave to get the time limits and use
% that range for all conditions. The limits are set when the wave crosses
% the threshold upwards and downwards

% 'auto_fixed75%' set the frame when the maximum condition reaches the 75%
% of the peak

% 'auto_fixed100%' set the frame when the maximum condition reaches the peak


% FOR WAVES
%---------------------------------------------------------------
% Note that for the computation initial and ending times from tiles and
% early peaks, the dF movie is always used
wave_units = 'dF'; % 'dF' '%F'
draw_thresh = 1;

% MANUAL PARAMETERS (if 'manual' options are selected)
%---------------------------------------------------------------

% if thresh_mode = 'manual' SET:
manual.clims = [];
manual.thresh = [];

% if timerange_mode = 'manual' SET:
manual.start_ms = [];
manual.end_ms = [];

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on=3;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 0; %
setting.forcein = 0; %

% NEWMAP 1: from jet
%     n = 512;
%     bigjet = jet(n);
%     cmap_tile = bigjet(round(n*0.3):end, :);


% NEWMAP 2: from new function

% cmap_tile = colormap_loadBV4(256);
cmap_tile = colormap_loadBV2(256);
%     cmap_tile = colormap_loadBV();

% SUBJECTS SELECTION (ROWSA FROM THE 'groupplot' CELL STRUCTURE
% FOR NO-PAIN
...................
    
% FOR PAIN
...................
    if strcmpi( analysisref, 'new4_group17_pulsito') % del panel original: para dm2+dm4 vs blank
    sel_rows = [24 25 26 29 30];
    elseif strcmpi( analysisref, 'new4_group19_dolor_n5') % del panel original: para dm2+dm4 vs blank
        sel_rows = [31:35];
    elseif strcmpi( analysisref, 'new4_group19_dolor_n4') % del panel original: para dm2+dm4 vs blank
        sel_rows = [31:34];
        
        % FOR INTENSITY EFFECT BOXPLOTS
        ...................
            % from: load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot5_boxplot_Ieffect.mat') % ATTENTION: select cases manually (there are repetitions)
        
    elseif strcmpi( analysisref, 'new4_group20_3intens') % del panel original: para dm2+dm4 vs blank (n=6)
        sel_rows = [11:16]; %n=4
        
        % FOR CHORRITO
        % ...................
        % from: load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot6_chorrito.mat') % ATTENTION: select cases manually (there are repetitions)
    elseif strcmpi( analysisref, 'new4_group16_chorrito_n3') % del panel original: para dm2+dm4 vs blank
        sel_rows = 1:3;
    elseif strcmpi( analysisref, 'new4_group16_chorrito_n4') % del panel original: para dm2+dm4 vs blank
        sel_rows = 1:4;
        
    elseif strcmpi( analysisref, 'TEST') % del panel original: para dm2+dm4 vs blank
        sel_rows = 31;
    end
    
    
    %% LOOP ACROSS FISH
    
    % preallocate output matrix

    for rowi  =  makeRow(sel_rows) %1:size(groupplot,1) SELECT included fish+condition
        
        
        % LOAD SETTINGS FROM STRUCTURE
        nfish = groupplot{rowi,1};
        cond_codes = groupplot{rowi,3};
        fact_thresh =  groupplot{rowi,4};
        
        %----------------------------------------------------------------
        % LOAD DATA
        %----------------------------------------------------------------
        VSDI = TORus ('load', nfish);
        VSDmov = TORus('loadmovie',nfish,movie_ref);
        
        %----------------------------------------------------------------
        % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
        %----------------------------------------------------------------
        
        rejectidx  = compute_rejectidx(VSDI, reject_on, setting);
        
        %----------------------------------------------------------------
        % SELECT ROI
        %----------------------------------------------------------------
        selroi =name2idx(selroinames, VSDI.roi.labels_circ);
        roilabels = VSDI.roi.labels_circ;
        masks =  VSDI.roi.circle.mask;
        
        %----------------------------------------------------------------
        % PEAK-FINDING FUNCTION PARAMETERS (to center the early-peak function)
        %----------------------------------------------------------------
        
        feedf.window.min = [-100 100]; %'feed-function' structure
        feedf.window.max = [0 1200]; % where to find the max peak
        feedf.window.movsum = 50; %ms
        feedf.window.basel = [-25 0]; %cambiar a  -100
        feedf.window.slope=50;
        feedf.window.wmean=[0 350];
        
        %----------------------------------------------------------------
        % trage for peak max-values finding
        %----------------------------------------------------------------
        
        trange = [0 1300]; % where to find the max values
        trange = dsearchn(VSDI.timebase,trange');
        
        trange = trange(1):trange(2);
        timebase_adj = VSDI.timebase(trange);
        
        %----------------------------------------------------------------
        % sample frequency for lowpass filter
        %----------------------------------------------------------------
        
        fs = 1000/(VSDI.info.stime);
        
        %% ----------------------------------------------------------------
        % CODE
        %----------------------------------------------------------------
        
        %................................................................
        % PRELIMINARY LOOP TO GET PARAMETERS
        %................................................................
        %----------------------------------------------------------------
        
        % fact_thresh = 0.4;
        
        % -------------------------------------------------------
        % GET REFERENCE-ROI MASK
        % -------------------------------------------------------
        
        idxDm4 =name2idx(refroiname, roilabels);
        ref_roimask = masks(:,:,idxDm4);
        
        %-------------------------------------------------------------------
        % GET MAX FROM ALL CONDITIONS TO SET AS REPRESENTATION'S THRESHOLD FOR
        % TILES
        %-------------------------------------------------------------------
        ci = 0;
        for condi = makeRow(cond_codes)
            ci = ci+1;
            
            % -------------------------------------------------------
            % SELECT TRIALS
            % -------------------------------------------------------
            
            [sel_trials] = find(VSDI.condition(:,1)==condi);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            % keep only first trials - for those cases in which only the first 'n'
            % trials are to be kept
            if keep_onlyNtrials
                sel_trials = sel_trials(1:keep_onlyNtrials);
            end
            
            back = VSDI.backgr(:,:,VSDI.nonanidx(1));
            %to plot single trial
            movie2plot = mean(VSDmov.data(:,:,trange,sel_trials),4);
            movie2plot(:,:,end+1) = back; %clean non-blured background
            
            % -------------------------------------------------------
            % GET MAX-VALUES FOR BOTH METHODS
            % -------------------------------------------------------
            % for 'moviebased_thresh_max'
            tempmax = movmean(movie2plot(:,:,1:end-1),5); %temporal smooth to get the max value par: 3
            
            allpix_localmax(ci) = max(tempmax(:));
            
            % For 'wavebased_thresh_local' and 'wavebased_thresh_max'
            dm4_wave = roi_TSave(movie2plot,ref_roimask);
            %     dm4_wave = movmean(dm4_wave ,5);
            dm4_wave = lowpass(dm4_wave,10,fs);
            
            dm4wave_localmax(ci) = max(dm4_wave);
            
            clear tempmax movie2plot  sel_trials dm4_wave
            
        end % for condi
        
        
        % CHECK POINT - whether the maximum calculated from waves is the same from
        % the movie
        
        if max(allpix_localmax) ~= max(dm4wave_localmax)
            warning('there seems to be a maximum in a region outside the roi of interest - that might yield and error or undesired results.')
        end
        
        % GET MAX VALUE FOR THRESHOLD CALCULATION LATER IN THE CODE
        switch thresh_mode
            case 'moviebased_thresh_max'
                [maxval, maxidx] = max(allpix_localmax); % the idx gives the condition with maximum movie-based value %BASED ON MOVIE
                
                warning('The MAX VALUE used to calculate the threshold is based on the maximum movie-pixel.')
                
            case {'wavebased_thresh_max' 'wavebased_thresh_local'}
                [maxval, maxidx] = max(dm4wave_localmax); % the idx gives the condition with maximum movie-based value %BASED ON MOVIE
                warning('The MAX VALUE used to calculate the threshold is based on the maximum reached by dm4 wave.')
        end
        
        
        % -------------------------------------------------------
        % FOR MODE '' GET dF WAVE OF MAXIMUM WAVE TO ESTABLISH A FIXED START-END TIME
        % -------------------------------------------------------
        
        ... GET AVERAGE MOVIE
            %----------------------------------------------------------------
        [sel_trials] = find(VSDI.condition(:,1)==cond_codes(maxidx));
        
        if reject_on
            sel_trials = setdiff(sel_trials, rejectidx);
        end
        
        
        % keep only first trials - for those cases in which only the first 'n'
        % trials are to be kept
        if keep_onlyNtrials
            sel_trials = sel_trials(1:keep_onlyNtrials);
        end
        
        
        back = VSDI.backgr(:,:,VSDI.nonanidx(1));
        
        %to plot single trial
        maxmovie = mean(VSDmov.data(:,:,trange,sel_trials),4);
        maxmovie(:,:,end+1) = back; %clean non-blured background
        
        % CALCULATE max WAVE FOR dm4 TO ESTABLISH START-END TIMES
        % -------------------------------------------------------
        max_wave = roi_TSave(maxmovie,ref_roimask);
        % max_wave = movmean(max_wave ,5);
        max_wave = lowpass(max_wave,10,fs);
        
        
        clear movie2plot  back sel_trials
        
        % ----------------------------------------------------------------
        % LOOP TO GET TILES, EARLY-PEAK FRAMES AND WAVES
        %----------------------------------------------------------------
        ci = 0;
        
        for condi = makeRow(cond_codes)
            ci = ci+1;
            %----------------------------------------------------------------
            ... GET AVERAGE MOVIE
                %----------------------------------------------------------------
            [sel_trials] = find(VSDI.condition(:,1)==condi);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            % keep only first trials - for those cases in which only the first 'n'
            % trials are to be kept
            if keep_onlyNtrials
                sel_trials = sel_trials(1:keep_onlyNtrials);
            end
            
            
            back = VSDI.backgr(:,:,VSDI.nonanidx(1));
            
            %to plot single trial
            movie2plot = mean(VSDmov.data(:,:,trange,sel_trials),4);
            movie2plot(:,:,end+1) = back; %clean non-blured background
            
            % -------------------------------------------------------
            % CALCULATE dF WAVE FOR dm4 TO ESTABLISH START-END TIMES
            % -------------------------------------------------------
            
            dm4_wave = roi_TSave(movie2plot,ref_roimask);
            %     dm4_wave = movmean(dm4_wave ,5);
            dm4_wave = lowpass(dm4_wave,10,fs);
            wave_timebase = timebase_adj;
            
            %----------------------------------------------------------------
            % TILES SETTINGS (for 'tiles' and early-peak'
            %----------------------------------------------------------------
            
            tileset.nrowcol = [1 6];
            tiles.backgr = back;
            tileset.interp =6;
            
            
            fr_pre = 3; %nº of frames before rising the threshold (to set the initial time)
            fr_post = 2;
            
            
            %----------------------------------------------------------------
            % GET REPRESENTATION THRESHOLD ACCORDING TO THE MODE
            %----------------------------------------------------------------
            
            switch thresh_mode
                case 'moviebased_thresh_max'
                    tileset.clims = [0 maxval*fact_clim];
                    tileset.thresh = [-maxval*fact_thresh maxval*fact_thresh];
                case 'wavebased_thresh_local'
                    localmax = dm4wave_localmax(ci);
                    maxwave = max(dm4wave_localmax);
                    tileset.clims = [0 maxwave*fact_clim]; % the colormap has to be the same  for all conditions to be comparable
                    tileset.thresh = [-localmax*fact_thresh localmax*fact_thresh];
                case 'wavebased_thresh_max'
                    maxwave = max(dm4wave_localmax);
                    tileset.clims = [0 maxwave*fact_clim]; % the colormap has to be the same for all conditions to be comparable
                    tileset.thresh = [-maxwave*fact_thresh maxwave*fact_thresh];
                    
                case 'manual'
                    tileset.clims = manual.clims;
                    tileset.thresh = manual.thresh;
                    fact_clim = []; %to avoid fake values when entering the title
                    fact_thresh = [];
            end
            
            % GET START/END/PEAK IDX (and dummy values in case that
            % the wave does not reach the threshold, that is if idx.start = [])
            %----------------------------------------------------------------
            switch timerange_mode
                case 'auto_localwave'
                    idx.start = find(dm4_wave > tileset.thresh(2), 1, 'first')- fr_pre;
                    idx.end= find(dm4_wave > tileset.thresh(2), 1, 'last') + fr_post;
                case 'auto_fixed'
                    % gets as timerange the time that the maximum reference wave
                    % among conditions takes to reach the threshold. Makes sense
                    % when comparing various conditions
                    idx.start = find(max_wave > tileset.thresh(2), 1, 'first')- fr_pre;
                    idx.end= find(max_wave > tileset.thresh(2), 1, 'last') + fr_post;
                case 'auto_fixed75%'
                    %                 temp = devo_peak2peak(max_wave, wave_timebase, feedf.window, [],  'movsum' , 0, 0);
                    %                 temppeak =temp.peakminusbasel;
                    % DEPRECATED: the max value has to be considered for the threshold (as done
                    % before). If this is consideredworng, then it has to be switched where
                    % 'tileset.thresh' is computed
                    temppeak = tileset.thresh(2)/fact_thresh;
                    idx.start = find(max_wave > tileset.thresh(2), 1, 'first');
                    idx.end= find(max_wave > temppeak*0.75, 1, 'first');
                    clear temppeak
                    
                case 'auto_fixed100%'
                    % the final time is the peak time of the larger wave among
                    % conditions
                    temppeak = tileset.thresh(2)/fact_thresh;
                    idx.start = find(max_wave > tileset.thresh(2), 1, 'first');
                    idx.end= find(round(max_wave,4) >= round(temppeak,4), 1, 'first');
                    clear temppeak
                    
                case 'manual'
                    idx.start = []; %it'll be assigned later in the code
                    idx.end = [];
            end
            
            if  ~ isempty (idx.start)
                
                temp = devo_peak2peak(dm4_wave, timebase_adj, feedf.window, [],  'movsum' , 0, 0);
                idx.peak =temp.peakidx(2);
                clear temp
                % warning ('THE CODE HAS BEEN HACKED BY CHANGING THE idx.peak - REVERSE WHEN YOU FINISH')
                %             idx.peak =temp.peakidx(2)-20;
                
            elseif isempty (idx.start)
                
                idx0 = dsearchn(timebase_adj, 0);
                ntiles = tileset.nrowcol(1) * tileset.nrowcol(2);
                
                idx.start = idx0;
                
                idx.end = idx0 +ntiles;
                idx.peak = idx0 +1;
            end
            
            
            %----------------------------------------------------------------
            % TILES
            %----------------------------------------------------------------
            
            if plottiles
                
                switch timerange_mode
                    case 'auto_localwave'
                        tileset.start_ms = timebase_adj(idx.start);
                        tileset.end_ms = timebase_adj(idx.end);
                        %         tileset.peak_ms = timebase_adj(idx.peak);
                    case 'manual'
                        tileset.start_ms = dsearchn(timebase_adj, manual.start_ms);
                        tileset.start_ms = (timebase_adj(tileset.start_ms));
                        
                        tileset.end_ms = dsearchn(timebase_adj, manual.end_ms);
                        tileset.end_ms = (timebase_adj(tileset.end_ms));
                        
                    case {'auto_fixed', 'auto_fixed75%' , 'auto_fixed100%'}
                        % gets as timerange the time that the maximum reference wave
                        % among conditions takes to reach the threshold. Makes sense
                        % when comparing various conditions
                        tileset.start_ms = timebase_adj(idx.start);
                        tileset.end_ms = timebase_adj(idx.end);
                        
                        
                end
                
                %         cmap_tile = colormap_loadBV2(256);
                imtiles = plot_tilemovie_custom_interp2(movie2plot, timebase_adj, tileset, [], cmap_tile);
                
                nframes = tileset.nrowcol(1)*tileset.nrowcol(2); %for the title
                
                titulo = [num2str(VSDI.ref) 'cond' num2str(condi) '.clim=' num2str(round(tileset.clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(tileset.thresh(2),1)) '(' num2str(fact_thresh*100) '%)'];
                sgtitle(titulo)
                
                savename= ['TILES_' num2str(VSDI.ref) movie_ref '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' thresh_mode  '_fact' num2str(fact_thresh) '_clim'  num2str(fact_clim) '_Trange_' timerange_mode 'Pre' num2str(fr_pre) 'Post' num2str(fr_post) 'Fr' num2str(nframes) movie_ref '.jpg'];
                
                if savetiles
                    %                         saveas(gcf, fullfile(savein, savename), 'jpg')
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                    print(fullfile(savein, savename),'-r600','-djpeg') % prints it as you see them
                    close
                end
                
                clear tileset.start_ms tileset.end_ms
                
            end
            
            % ----------------------------------------------------------------
            % PIXELS ABOVE THRESHOLD (for the last tiles frame)
            % ----------------------------------------------------------------
            
            % Get logic map of pixels above threshold
            abovethresh = imtiles(:,:,end) > tileset.thresh(2);
            
            % Intersect each roi with the roi mask
            for n = 1:numel(selroi)
                
                % Get anatomical- roi names from circle-roi
                Aroiname  = selroinames{n};
                Aroiname = Aroiname(1:end-1);
                % load anatomical masks
                Aidx = find(strcmpi(Aroiname, VSDI.roi.labels));
                maskA = VSDI.roi.manual_mask(:,:,Aidx);
                
                % get intersection
                pixU = abovethresh & maskA ;
                
                % STORE TO LATER PRINT IN OUTPUT EXCEL
                % .............................
                
                abovepix(rowi,n) = sum(pixU(:));
                
                npixroi = sum(maskA(:));
                abovepixP(rowi,n) = round((sum(pixU(:))/npixroi)*100);
                
                
                
                % PLOT CONTROL FIGURE
                % .............................
                f = figure;
                
                subplot(1,3,1)
                image(imfuse(VSDI.backgr(:,:,VSDI.nonanidx(1)), abovethresh))
                title('pix above thresh')
                axis image
                
                subplot(1,3,2)
                image(imfuse(VSDI.backgr(:,:,VSDI.nonanidx(1)), maskA))
                title('roi mask')
                axis image
                
                subplot(1,3,3)
                image(imfuse(maskA, abovethresh))
                axis image
                title(['inters:' num2str(sum(pixU(:))) 'pix (' num2str(abovepixP(rowi,n)) '%)'])
                tit = ['INTERSECT' num2str(VSDI.ref) '-' num2str(condi) '- fact' num2str(tileset.thresh(2)) '-' Aroiname];
                saveas(f, fullfile(savein,[tit '.jpg']), 'jpg')
                close
                
            end
            
            % Add number of fish and condition
            abovepix(rowi,n+1) = VSDI.ref;
            abovepix(rowi,n+2) = condi;
            abovepix(rowi,n+3) = fact_thresh;

            % Add number of fish and condition
            abovepixP(rowi,n+1) = VSDI.ref;
            abovepixP(rowi,n+2) = condi;
            abovepixP(rowi,n+3) = fact_thresh;
            
            %----------------------------------------------------------------
            ... PLOT EARLY-PEAK
                %----------------------------------------------------------------
            
            % GET REPRESENTATION THRESHOLD ACCORDING TO THE MODE
            %----------------------------------------------------------------
            
            %     if isempty (tileset.start_ms)
            %         error('The threshold is too high (lower the value of "fact_thresh"')
            %     end
            
            
            if plot_earlypeak
                tileset.nrowcol = [1 2];
                
                tileset.start_ms = timebase_adj(idx.start);
                tileset.end_ms = timebase_adj(idx.peak);
                
                % HACKING SNIPPET !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % %         to plot specific frame instead of real peak
                %         warning ('THE CODE HAS BEEN HACKED BY CHOOSING A SPECIFIC FRAME - REVERSE WHEN YOU FINISH')
                %         tileset.end_ms = dsearchn(timebase_adj, 550);
                %         tileset.end_ms = (timebase_adj(tileset.end_ms));
                % END OF HACKING SNIPPET !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                nframes = 2;
                plot_tilemovie_custom_interp2(movie2plot, timebase_adj, tileset, [], cmap_tile);
                
                titulo = [num2str(VSDI.ref) '_cond' num2str(condi) '.clim=' num2str(round(tileset.clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(tileset.thresh(2),1)) '(' num2str(fact_thresh*100) '%)'];
                sgtitle(titulo)
                
                savename= ['EARLYnPEAK_frames' num2str(VSDI.ref) '_' thresh_mode  'trials_factors' num2str(fact_thresh)  '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) '-' num2str(fact_clim) 'pre' num2str(fr_pre) 'post' num2str(fr_post) 'fr' num2str(nframes) movie_ref '.jpg'];
                
                if save_earlypeak
                    %                         saveas(gcf, fullfile(savein, savename), 'jpg')
                    set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                    print(fullfile(savein, savename),'-r600','-djpeg') % prints it as you see them
                    close
                end
                clear tileset.start_ms tileset.end_ms
                
            end
            
            %     %----------------------------------------------------------------
            %     % CONTOURS - STORE TO LATER PLOT
            %     %----------------------------------------------------------------
            %     if plot_contour
            %
            %         backgr = movie2plot(:,:,end);
            %         if highdefinition
            %             backgr = interp2(backgr,5, 'cubic');
            %         end
            %
            %         % GET FRAME TO INTERSECT
            %         switch contour_of
            %
            %             case 'peak'
            %                 % .....................................................................
            %                 ... PEAK CONTOUR (depending on 'thresh_mode')
            %                 % .....................................................................
            %                 tileset.peak = timebase_adj(idx.peak);
            %                 if isempty (tileset.peak)
            %                     error('The threshold is too high to get all contours (lower the value of "fact_thresh"')
            %                 end
            %
            %                 frame = movie2plot(:,:,idx.peak);
            %
            %             case 'early'
            %                 % .....................................................................
            %                     ... EARLY CONTOUR (depending on 'thresh_mode' & 'timerange_mode' )
            %                 % .....................................................................
            %                 if isempty (idx.start)
            %                     error('The threshold is too high to get all contours (lower the value of "fact_thresh"')
            %                 end
            %
            %                 frame = movie2plot(:,:,idx.start);
            %
            %
            %         end
            %
            %                 if highdefinition
            %                     frame = interp2(frame,5, 'cubic');
            %                 end
            %
            %                 % TURN FRAME INTO SURFACE
            %                 for x = 1:size(frame,1)
            %                     for y = 1:size(frame,2)
            %                         xs(x,y) = x;
            %                         ys(x,y) = y;
            %                         zs(x,y) = frame(x,y);
            %                     end
            %                 end
            %
            %                 % INTERSECT WITH PLANE
            %                 thresh  =tileset.thresh;
            %                 figure
            %                 [M,~] = contourf( xs,ys,zs, [thresh thresh]);
            %                 close;
            %                 %TURN INTO COORDS
            %                 if ~isempty(M)
            %                     [x,y,z] = C2xyz(M) ; %get coordinates of the intersection
            %                     structC.n= length(x);
            %                     structC.X = x;
            %                     structC.Y = y;
            %                     structC.Z = z;
            %
            %                 else
            %                     structC.n= 0;
            %                     structC.X = 0;
            %                     structC.Y = 0;
            %                     structC.Z = 0;
            %                 end
            %                 clear x y z M
            %
            %                 contours{ci} = structC;
            %
            %     end % if plot_contour
            
            clear idx
            
            
            
            
            
        end % for condi
        % blob()
        
        %% ----------------------------------------------------------------
        % EXPORT TO EXCEL
        %----------------------------------------------------------------
        
        % Create parameters 
        params{1,1} = thresh_mode;
        params{2,1} = timerange_mode;
        params{3,1} = 'refroi' ; params{3,2}  =refroiname;
        params{4,1}= analysisref;
        params{6,1} = date;
        params{6,3} = ['source:'  mfilename('fullpath')];


        % Wrte each fish - condition in a separate sheet
        labels = VSDI.roi.labels_circ(selroi);
        labels{end + 1} = 'fish';
        labels{end + 1} = 'cond';
        labels{end + 1} = 'thresh';

            excelname = fullfile(savein, [analysisref 'PIXELS' movie_ref '_rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls']);
            
            for col = 1:numel(labels)
                labels{1,col}= labels{col};
            end
            
            % write output (new sheet for each fish
            writematrix (abovepix, excelname, 'sheet',  'npix')
            writematrix (abovepixP, excelname, 'sheet',  '%pix')
            
            writecell (labels, excelname, 'sheet', 'labels')
            writecell (params, excelname, 'sheet', 'params')

            %writecell (params, excelname, 'sheet','param')
      %%      
    end
    
    %
    % %% ----------------------------------------------------------------
    % % PLOT AND SAVE CONTOURS
    % %----------------------------------------------------------------
    % if plot_contour
    %
    %     timecmap = lines(numel(cond_codes));
    %
    %     figure
    %     imagesc(backgr); colormap bone
    %     axis image
    %     hold on
    %
    %     ci = 0;
    %     for condi = makeRow(cond_codes)
    %         ci = ci+1;
    %
    %         cont = contours{ci};
    %
    %         if cont.n > 0
    %             for n= 1:cont.n
    %                 x = cont.X{n};
    %                 y= cont.Y{n};
    %
    %                 plot(y, x, 'linewidth' , 1.5, 'color', timecmap(ci,:), 'displayname', num2str(condi)); hold on
    %                 clear x y
    %
    %             end
    %         end
    %
    %     end
    %
    %     %         sgtitle([num2str(VSDI.ref), movie_ref, 'rej', num2str(reject_on),'-', VSDI.info.Sside '(threh' num2str(fact_thresh) '%max)'])
    %     %         savename= ['map_peak_' num2str(VSDI.ref) movie_ref  VSDI.info.Sside '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials'];
    %
    %     legend ('Location','northeastoutside')
    %     sgtitle([num2str(VSDI.ref), 'rej', num2str(reject_on),'-', thresh_mode '- (threh' num2str(fact_thresh) '%max)' 'trange:' timerange_mode])
    %
    %     savename= ['CONTOUR' contour_of '_' num2str(VSDI.ref) movie_ref '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' thresh_mode '_fact' num2str(fact_thresh)];
    %
    % end
    %
    % if save_contour
    %     saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
    %
    % %                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
    % %                 print(fullfile(savein, [savename '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them %STILL TO TEST!
    %
    %     close
    % end
    %
    % % end % FOR THE 'rowi' LOOP - silence when the loop is not used
    %
    % blob()
    % % newsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes_code/03_figure_sketch/def_figs/tiles';
    % % saveas(gcf, fullfile(newsave,savename) 'jpg')
    %
    % %----------------------------------------------------------------
    % % PRINT INCLUDED AND EXCLUDED TRIALS
    % % ----------------------------------------------------------------
    % %
    % % for condi = cond_codes
    % %
    % %     [sel_trials] = find(VSDI.condition(:,1)==condi);
    % %     local_reject = intersect(sel_trials, rejectidx);%just to display later;
    % %     sel_trials = setdiff(sel_trials, rejectidx);
    % %
    % %     disp(['Included trials for condition' num2str(condi) ':' ])
    % %     disp(sel_trials)
    % %     disp('%')
    % %     disp(['Rejected trials for condition' num2str(condi) ':'])
    % %     disp(local_reject)
    % % end
    
    % Updates
    % FORK to get a measure of pixels above threshold
    % 21/10/22_ add 'keep...
    % 17/10/22:
    % -add 'autofixed100%'
    % -remove 'fr_pre' substraction from cases.
    % -Fixed bug related to the maxval and add the code that calculates it appropiately for each 'thresh_mode'
    % 'auto_fixed75%' and 'auto_fixed100%'
    % - use 'trange' and 'timebase_adj' instead of '1:end' and 'VSDI.timebase'
    % 26/09/22: Fix bug: [sel_trials] = find(VSDI.condition(:,1)==cond_codes(maxcondi));
    % 24/09/22: Update plot early peak
    % 21/09/22: add timerange_mode to the code
    % 18/09/22: adapted from:
    % '/home/tamara/Documents/MATLAB/VSDI/MOT1x/plot_code/tiles_waves_earlypeak_contour.m'
    % 17/09/22: add 'manual' option for thresh_mode
    % 12/03/22 : compute_rejectidx function and improve 'cmap'