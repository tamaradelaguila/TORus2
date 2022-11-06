%% CODE: TILES + WAVEs (all roi for each condition)

clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

%----------------------------------------------------------------
% @SET: fish + configuration
%----------------------------------------------------------------
%
% nfish = 18;
% cond_codes = [501]; % DO NOT INCLUDE THE BLANK CONDITION (or the threshold will fail)
% ATT: the threshold will be computed respect to the maxim um condition


% for rowi = [1:2]
% clearvars -except grouptiles rowi path
% TO LOOP FROM STRUCTURE
...................................
    % % load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/grouptiles_all.mat')
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/grouptiles_all.mat')
%
% nfish = grouptiles{rowi,1};
% cond_codes = grouptiles{rowi,2};
% selroinames = grouptiles{rowi,3};
%  ref_croproi = grouptiles{rowi,4};
% ....................................
 

% TO SKIP THE LOOP
% ....................................
nfish = 11;
cond_codes = [401 402 404];
% selroinames = {'dm4m_R2',  'dm2_R2' ,'dm1_R','dldm_R'};%dm3 ORIGINAL

ref_croproi = 'dm4_R'; % respect which the start/end time of the tiles will be set
% ....................................


plot_contour = 1;
save_contour = 1;
highdefinition = 0;

contour_of = 'early'; % 'early' 'peak'
% Both codes find the thresh/frame that meet the condition of npix above
% thresh:
% ...early - threshold is fixed, finds the frame 
% ...peak - frametime is fixed at peak level, finds the thresh 
% condition

plottiles = 1; %also plots early-peak
savetiles = 1;


movie_ref = '_18filt6'; % input movie '_17filt5'
% movie_ref = '_17filt5'; % input movie '_17filt5'
% movie_ref = '_15filt5'; % usado para 220608


% savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/def_figs/tiles/'; %
savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/'; %

% FOR TILES AND CONTOURS
%---------------------------------------------------------------
fact_thresh = 0.35; % @SET : limits parameters
fact_clim= 1.5;

thresh_mode = 'moviebased_thresh_local'; % 'moviebased_thresh_max' 'moviebased_thresh_local'.
% timerange_mode = 'auto_localwave' ; %'auto_localwave', 'manual', 'auto_fixed'

% for contours: guiding parameters 
npix = 10; %set how many pix above thresh
nfr = 5; %for how many nconsecutive frames
startfrom_ms = 0; %which index use as first
% auto_localwave - sets the timerange according to the condition's
% reference wave

% auto_fixed - when there several conditions, it calculates the reference wave of
% all conditions, takes the highest wave to get the time limits and use that range for all conditions


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
reject_on= 0;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.forcein = 0; %


% NEWMAP 2: from new function

cmap_tile = colormap_loadBV2(256);

%% LOAD AND SETTINGS
%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus ('load', nfish);
VSDmov = TORus('loadmovie',nfish,movie_ref);

% SELECT RECTANGULAR ROI TO CROP MOVIE
%----------------------------------------------------------------
roilabels = VSDI.roi.labels_rect;

selroi =name2idx( ref_croproi,roilabels);
masks = VSDI.roi.rect.mask;

%----------------------------------------------------------------
% CROP MOVIES
%----------------------------------------------------------------
cropmask = masks(:,:,selroi);

% cropvoxels = repmat(cropmask,  1,1, nframe); % 3D mask to select only the pixels from the mask
% d = size(cropvoxels);

for tri = 1: size(VSDmov.data, 4)
    for fri = 1:size(VSDmov.data,3)
        tempframe= VSDmov.data(:,:,fri,tri);
        cropmovie(:,:,fri,tri) =tempframe.*cropmask;
        clear tempframe
    end
    cropbackgr(:,:,tri) = VSDI.backgr(:,:,tri).*cropmask;
    clear tempmov
end

%----------------------------------------------------------------
% THIS DOES NOT WORK
% for tri = 1:ntri
%     tempmov = VSDmov.data(:,:,:,ntri);
%     for fram = 1:nframe
% cropmovie2(:,:,:,tri) =tempmov(cropvoxels);
%     end
% clear tempmov
% end


%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------

rejectidx  = compute_rejectidx(VSDI, reject_on, setting);

%----------------------------------------------------------------
% PEAK-FINDING FUNCTION PARAMETERS (to center the early-peak function)
%----------------------------------------------------------------

feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 1200]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-25 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350];

feedf.method = 'movsum';
%% ----------------------------------------------------------------
% CODE
%----------------------------------------------------------------
%----------------------------------------------------------------
% PRELIMINARY LOOP TO GET PARAMETERS
%----------------------------------------------------------------
%----------------------------------------------------------------

% fact_thresh = 0.4;

% -------------------------------------------------------
% GET REFERENCE-ROI MASK ----UNUSED
% -------------------------------------------------------

% idxDm4 =name2idx( ref_croproi, roilabels);
% ref_roimask = masks(:,:,idxDm4);

%%
%-------------------------------------------------------------------
% GET MAX FROM ALL CONDITIONS TO SET AS REPRESENTATION'S THRESHOLD FOR
% TILES
%-------------------------------------------------------------------
switch thresh_mode
    case 'moviebased_thresh_max'
        ci = 0;
        for condi = makeRow(cond_codes)
            ci = ci+1;
            
            [sel_trials] = find(VSDI.condition(:,1)==condi);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            
            back = cropbackgr(:,:,VSDI.nonanidx(1));
            
            
            %to plot single trial
            movie2plot = mean(cropmovie(:,:,:,sel_trials),4);
            movie2plot(:,:,end) = back; %clean non-blured background
            
            % -------------------------------------------------------
            % GET MAX-VALUES FOR BOTH METHODS
            % -------------------------------------------------------
            % for 'moviebased_thresh_max'
            tempmax = movmean(movie2plot(:,:,1:end-1),5); %temporal smooth to get the max value par: 3
            
            maxall(ci) = max(tempmax(:));
            
            
        end
        
        [maxval, maxidx] = max(maxall); % the idx gives the condition with maximum value
        
        % -------------------------------------------------------
        % FOR MODE '' GET dF WAVE OF MAXIMUM WAVE TO ESTABLISH A FIXED START-END TIME
        % -------------------------------------------------------
        maxcondi = maxidx;
        
        ... GET AVERAGE MOVIE
            %----------------------------------------------------------------
        [sel_trials] = find(VSDI.condition(:,1)==cond_codes(maxcondi));
        
        if reject_on
            sel_trials = setdiff(sel_trials, rejectidx);
        end
        
        back = cropbackgr(:,:,VSDI.nonanidx(1));
        
        %to plot single trial
        maxmovie = mean(cropmovie(:,:,:,sel_trials),4);
        maxmovie(:,:,end) = back; %clean non-blured background
        
        clear movie2plot  back sel_trials
        
end % switch

%%
% ----------------------------------------------------------------
% LOOP TO GET TILES, EARLY-PEAK FRAMES AND WAVES
%----------------------------------------------------------------
ci = 0;

for condi = makeRow(cond_codes)
    ci = ci+1;
    %----------------------------------------------------------------
    ... GET AVERAGE MOVIE
        %----------------------------------------------------------------
    [sel_trials] = find(VSDI.condition(:,1)==condi)';
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    
    back = cropbackgr(:,:,VSDI.nonanidx(1));
    
    %to plot single trial
    movie2plot = mean(cropmovie(:,:,:,sel_trials),4);
    movie2plot(:,:,end) = back; %clean non-blured background
    
    %----------------------------------------------------------------
    % GET REPRESENTATION THRESHOLD ACCORDING TO THE MODE
    %----------------------------------------------------------------
    
    switch thresh_mode
        case 'moviebased_thresh_max'
            clims = [0 maxval*fact_clim];
            thresh =  maxval*fact_thresh;
            
        case 'manual'
            clims = manual.clims;
            thresh = manual.thresh;
            fact_clim = []; %to avoid fake values to enter the title
            fact_thresh = [];
            
        case 'moviebased_thresh_local'
            
            
            % for 'moviebased_thresh_max'
            localmov = movmean(movie2plot(:,:,1:end-1),5); %temporal smooth to get the max value par: 3
            
            [localmax, maxidx] = max(localmov(:));
            
            sz = size(localmov);
            ii = 0;
            for x = 1:sz(1)
                for y = 1:sz(2)
                    ii = ii+1;
                    pixwave = localmov(x,y,:);
                   temp = devo_peak2peak(pixwave, VSDI.timebase, feedf.window, [],feedf.method, 0, 0);
                    allpeak(ii) = temp.peakminusbasel;
                    allidx(ii) = temp.peakidx(2); 
                   clear pixwave temp
                end
            end
            clear ii
            peakval = max(allpeak);
            peakidx = max(allidx);
            
            % GET CLIMS AND THRESHOLD 
            
            clims = [0 localmax*fact_clim];
            thresh = localmax*fact_thresh;
            
            clear localmov localmax
            
        case 'moviebased_thresh_local'
            %base the threshold in wave

    end
    
    
    %     %----------------------------------------------------------------
    %     % TILES
    %     %----------------------------------------------------------------
    %
    %     if plottiles
    %
    %         switch timerange_mode
    %             case 'auto_localwave'
    %                 tileset.start_ms = VSDI.timebase(idx.start);
    %                 tileset.end_ms = VSDI.timebase(idx.end);
    %                 %         tileset.peak_ms = VSDI.timebase(idx.peak);
    %             case 'manual'
    %                 tileset.start_ms = dsearchn(VSDI.timebase, manual.start_ms);
    %                 tileset.start_ms = (VSDI.timebase(tileset.start_ms));
    %
    %                 tileset.end_ms = dsearchn(VSDI.timebase, manual.end_ms);
    %                 tileset.end_ms = (VSDI.timebase(tileset.end_ms));
    %
    %             case 'auto_fixedwave'
    %                 % gets as timerange the time that the maximum reference wave
    %                 % among conditions takes to reach the threshold. Makes sense
    %                 % when comparing various conditions
    %                 tileset.start_ms = VSDI.timebase(idx.start);
    %                 tileset.end_ms = VSDI.timebase(idx.end);
    %         end
    %
    % %         cmap_tile = colormap_loadBV2(256);
    %         plot_tilemovie_custom_interp2(movie2plot, VSDI.timebase, tileset, [], cmap_tile);
    %
    %         nframes = tileset.nrowcol(1)*tileset.nrowcol(2); %for the title
    %
    %         titulo = [num2str(VSDI.ref) 'cond' num2str(condi) '.clim=' num2str(round(clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(thresh(2),1)) '(' num2str(fact_thresh*100) '%)'];
    %         sgtitle(titulo)
    %
    %             savename= ['TILES_' num2str(VSDI.ref) movie_ref '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' thresh_mode  '_fact' num2str(fact_thresh) '_clim'  num2str(fact_clim) '_Trange_' timerange_mode 'Pre' num2str(fr_pre) 'Post' num2str(fr_post) 'Fr' num2str(nframes) movie_ref '.jpg'];
    %
    %         if savetiles
    %             %                         saveas(gcf, fullfile(savein, savename), 'jpg')
    %             set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
    %             print(fullfile(savein, savename),'-r600','-djpeg') % prints it as you see them
    %             close
    %         end
    %
    %         clear tileset.start_ms tileset.end_ms
    %
    %     end
    %
    
    
    %----------------------------------------------------------------
    % CONTOURS - STORE TO LATER PLOT
    %----------------------------------------------------------------
    if plot_contour
        
        if highdefinition
            back = interp2(back,5, 'cubic');
        end
        
        %----------------------------------------------------------------
        % GET FRAME TO INTERSECT
        %----------------------------------------------------------------
        
        switch contour_of
            
            case 'peak'
                % frame at peak level, the threshold is lowered until the
                % npix is reached
                
            peakthresh = peakval;
            step = peakthresh/100; 
            frameidx= peakidx;
                nabove = 0;
                stopflag = 0;
            
                while ~stopflag
                    frame = movie2plot(:,:,frameidx);
                    nabove = sum(frame(:)>= peakthresh);
                    
                    if nabove >= npix
                        stopflag = 1;
                    end
              peakthresh = peakthresh - step; 
                
                end

            thresh = peakthresh + step; %to compensate for the last addition
            realabove = nabove; % the name 'real' is inherited from the case 'early'
            
            fact_thresh = round((thresh)/peakval , 2); 
            
            case 'early'
                % Given a threshold, advances frame by frame until the npix
                % activation is met (also has checks that the activation is
                % mantained for a minimun number of frames)
                
                frameidx = dsearchn(VSDI.timebase, startfrom_ms) -1 ; %'-1' to compensate for the first addition in loop
                nabove = 0;
                stopflag = 0;
                %                 while nabove < npix
                %                     frameidx = frameidx+1;
                %                     frame = cropmovie(:,:,frameidx);
                %                     nabove = sum(frame(:)>= thresh);
                %
                %                 end
                %
                % % ............................................................
                % %  ONLY COUNTS PIXELS ABOVE THE THRESH AT ONCE
                % % ............................................................
                %                 while ~stopflag
                %                     frameidx = frameidx+1;
                %                     frame = movie2plot(:,:,frameidx);
                %                     nabove = sum(frame(:)>= thresh);
                %
                %                     for i= 1:nfr
                %                         tempframe = movie2plot(:,:,frameidx+i);
                %                         nabove = sum(tempframe(:)>= thresh);
                %
                %                         condition(i) = nabove> npix;
                %                     end
                %                     if sum(condition) == nfr
                %                        stopflag = 1;
                %
                %                     end
                %
                %                 end
                % ............................................................
                % CODE THAT ONLY COUNTS PIXELS ABOVE THE THRESH FOR A CERTAIN AMOUNT OF
                % NFRAMES
                % ............................................................

                while ~stopflag
                    frameidx = frameidx+1;
                    frame = movie2plot(:,:,frameidx);
                    nabove = sum(frame(:)>= thresh);
                    
                    % check which are really considered above thresh (if it
                    % stay for a while)
                    %                     realabove = 0;
                    
                    if nabove>npix
                        realabove = 0;
                        whichpix = find(frame(:)>= thresh);
                        sz= size(frame);
                        for i= 1:numel(whichpix)
                            pix = whichpix(i);
                            [row,col] = ind2sub(sz,pix); % GET COORDINATES OF THE PIXEL
                            wavesnip = squeeze(movie2plot(row,col,frameidx:frameidx+nfr-1));
                            if sum(wavesnip > thresh)==nfr
                                realabove = realabove +1;
                            end
                            %get wave, check if it meets the continuity and sum the pixel if so, sum to realabove
                        end
                        
                        if realabove> npix
                            stopflag = 1;
                            disp(['stop in frame' ,num2str(frameidx), '-', num2str(VSDI.timebase(frameidx)), 'ms'])
                        end
                        clear whichpix pix row col wavesnip
                    end
                    
                    % if the condition has not meet, assign a dummy value
                    % to 'frameidx'
                    if frameidx>= size(movie2plot,3)
                        frameidx = dsearchn(VSDI.timebase, startfrom_ms) -1 ; %'-1' to compensate for the first addition in loop
                        warning (['condition' num2str(condi) 'does not yield a result with the current parameters'])
                            break
                    end
                    
                end
        end
        
        if plottiles
            tileset.nrowcol= [1, 7];
            tileset.start_ms = VSDI.timebase(frameidx-1);
            tileset.end_ms = VSDI.timebase(frameidx+5);
            tileset.thresh = [-thresh thresh];
            tileset.clims = clims;
%             cmap_tile = colormap_loadBV2(256);
            
            plot_tilemovie_custom_interp2(movie2plot, VSDI.timebase, tileset, [], cmap_tile);
            
            
            titulo = [num2str(VSDI.ref) 'cond' num2str(condi) '.clim=' num2str(round(clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(thresh,1)) '(' num2str(fact_thresh*100) '%)' num2str(realabove) 'pixabove'];
            sgtitle(titulo)
            savename= ['TILES_' num2str(VSDI.ref) movie_ref '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' thresh_mode  '_fact' num2str(fact_thresh) '_clim'  num2str(fact_clim) num2str(realabove) 'pixabove_for' num2str(nfr) movie_ref '.jpg'];
            
            if savetiles
                  saveas(gcf, fullfile(savein, savename), 'jpg')
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%                 print(fullfile(savein, savename),'-r600','-djpeg') % prints it as you see them
                close
            end
            
            clear tileset.start_ms tileset.end_ms
            
            
        end %  if plottiles
        
        frame = movie2plot(:,:,frameidx);
        if highdefinition
            frame = interp2(frame,5, 'cubic');
        end
        
        % TURN FRAME INTO SURFACE
        for x = 1:size(frame,1)
            for y = 1:size(frame,2)
                xs(x,y) = x;
                ys(x,y) = y;
                zs(x,y) = frame(x,y);
            end
        end
        
        % INTERSECT WITH PLANE
        thresh  =thresh;
        figure
        [M,~] = contourf( xs,ys,zs, [thresh thresh]);
        close;
        %TURN INTO COORDS
        if ~isempty(M)
            [x,y,z] = C2xyz(M) ; %get coordinates of the intersection
            structC.n= length(x);
            structC.X = x;
            structC.Y = y;
            structC.Z = z;
            
        else
            structC.n= 0;
            structC.X = 0;
            structC.Y = 0;
            structC.Z = 0;
        end
        clear x y z M
        
        contours{ci} = structC;
        
    end % if plot_contour
    
    clear idx movie2plot rejectidx
    
end % for condi
% blob()


%% ----------------------------------------------------------------
% PLOT AND SAVE CONTOURS
%----------------------------------------------------------------
if plot_contour
    
    timecmap = lines(numel(cond_codes));
    
    figure
    imagesc(back); colormap bone
    axis image
    hold on
    
    ci = 0;
    for condi = makeRow(cond_codes)
        ci = ci+1;
        
        cont = contours{ci};
        
        if cont.n > 0
            for n= 1:cont.n
                x = cont.X{n};
                y= cont.Y{n};
                
                plot(y, x, 'linewidth' , 1.5, 'color', timecmap(ci,:), 'displayname', num2str(condi)); hold on
                clear x y
                
            end
        end
        
    end
    
    %         sgtitle([num2str(VSDI.ref), movie_ref, 'rej', num2str(reject_on),'-', VSDI.info.Sside '(threh' num2str(fact_thresh) '%max)'])
    %         savename= ['map_peak_' num2str(VSDI.ref) movie_ref  VSDI.info.Sside '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials'];
    
    legend ('Location','northeastoutside')
    sgtitle([num2str(VSDI.ref), 'rej', num2str(reject_on),'-', thresh_mode '- (threh' num2str(fact_thresh) '%max)' ])

  
    savename= ['CONTOUR' contour_of '_' num2str(VSDI.ref) movie_ref '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' thresh_mode '_fact' num2str(fact_thresh) '_' num2str(npix) 'pix-for' num2str(nfr) 'fr'];
    
end

if save_contour
    saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
    
    %                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
    %                 print(fullfile(savein, [savename '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them %STILL TO TEST!
    
    close
end

% end % FOR THE 'rowi' LOOP - silence when the loop is not used

blob()
% newsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes_code/03_figure_sketch/def_figs/tiles';
% saveas(gcf, fullfile(newsave,savename) 'jpg')

%----------------------------------------------------------------
% PRINT INCLUDED AND EXCLUDED TRIALS
% ----------------------------------------------------------------
%
% for condi = cond_codes
%
%     [sel_trials] = find(VSDI.condition(:,1)==condi);
%     local_reject = intersect(sel_trials, rejectidx);%just to display later;
%     sel_trials = setdiff(sel_trials, rejectidx);
%
%     disp(['Included trials for condition' num2str(condi) ':' ])
%     disp(sel_trials)
%     disp('%')
%     disp(['Rejected trials for condition' num2str(condi) ':'])
%     disp(local_reject)
% end

% Updates
% 26/09/22: Add 'npix' and 'nfr' conditions
% 24/09/22: Update plot early peak
% 21/09/22: add timerange_mode to the code
% 18/09/22: adapted from:
% '/home/tamara/Documents/MATLAB/VSDI/MOT1x/plot_code/tiles_waves_earlypeak_contour.m'
% 17/09/22: add 'manual' option for thresh_mode
% 12/03/22 : compute_rejectidx function and improve 'cmap'