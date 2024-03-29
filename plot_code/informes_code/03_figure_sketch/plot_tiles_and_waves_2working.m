close all
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

%----------------------------------------------------------------
% @SET: fish + conditions
%---------------------------------------------------------------
nfish = 11;
cond_codes = [401]; % DO NOT INCLUDE THE BLANK CONDITION (or the threshold will fail)
% ATT: the threshold will be computed respect to the maxim um condition

plottiles = 1; %also plots early-peak
savetiles = 0;

plot_earlypeak = 0;
save_earlypeak =0;

plotwaves=1;
savewaves = 0;

movie_ref = '_18filt6'; % input movie '_17filt5'
% movie_ref = '_17filt5'; % input movie '_17filt5'

% savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/def_figs/tiles/'; % 
savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/04_figure_sketch2/'; % 


% FOR WAVES
%---------------------------------------------------------------
% Note that for the computation initial and ending times from tiles and
% early peaks, the dF movie is always used 

waves_units = 'dF'; % 'dF' '%F'

% FOR TILES AND EARLY-PEAK FRAMES
%---------------------------------------------------------------

fact_thresh = 0.6; % @SET : limits parameters
fact_clim= 1.65;

thresh_mode = 'wavebased_thresh_max'; % 'moviebased_thresh_max' 'wavebased_thresh_max' 'wavebased_thresh_local' 'manual'. 'moviebased_thresh_max' is the one we have been using for tiles 
timerange_mode = 'refwave' ; %'refwave', 'manual' 


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

%----------------------------------------------------------------
% @SET: for roi-waves plot
%----------------------------------------------------------------

% selroinames = {'dm4m_R2',  'dm2_R2' ,'dm1_R','dldm_R2', 'dldr_R'};%dm3 ORIGINAL
% selroinames = {'dm4m_R2',  'dm2_R2' ,'dm1_R','dldm_R2'};%dm3 ORIGINAL

% selroinames = {'dm4m_R2',  'dm2_R2' ,'dm1_R', 'dldm_R2'};%dm3 % FROM '03figure_sketch'

selroinames = {'dm4m_R2',  'dm2_R2' ,'dm1_R', 'dldm_R'};%dm3

% selroinames = {'dm4m_R2', 'dm4m_L2',  'dm2_R2' , 'dm2_L2'};% IPSI VS CONTRA

refroiname = 'dm4m_R2'; % respect which the start/end time of the tiles will be set

roikind = 'circle'; %
% roikind = 'anat';

% NEWMAP 1: from jet
%     n = 512; 
%     bigjet = jet(n);
%     cmap_tile = bigjet(round(n*0.3):end, :);


% NEWMAP 2: from new function

    cmap_tile = colormap_loadBV2(256);
%     cmap_tile = colormap_loadBV();

%% LOAD / COMPUTE SETTINGS
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
if plotwaves
    switch roikind
        case 'circle'
            selroi =name2idx(selroinames, VSDI.roi.labels_circ);
            roilabels = VSDI.roi.labels_circ;
            masks =  VSDI.roi.circle.mask;
            
        case 'anat'
            selroi =name2idx(selroinames, VSDI.roi.labels);
            roilabels = VSDI.roi.labels;
            masks = VSDI.roi.manual_mask;
    end
    
else
    switch roikind
        case 'circle'
            roilabels = VSDI.roi.labels_circ;
            masks =  VSDI.roi.circle.mask;
            
        case 'anat'
            roilabels = VSDI.roi.labels;
            masks = VSDI.roi.manual_mask;
    end
    
end
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
%% CODE: TILES + WAVEs (all roi for each condition)
%----------------------------------------------------------------
    
% -------------------------------------------------------
    % GET REFERENCE-ROI MASK
    % -------------------------------------------------------
    
    idxDm4 =name2idx(refroiname, roilabels);
    roimask = masks(:,:,idxDm4);

%-------------------------------------------------------------------
% GET MAX FROM ALL CONDITIONS TO SET AS REPRESENTATION'S THRESHOLD FOR
% TILES
%-------------------------------------------------------------------
ci = 0;
for condi = makeRow(cond_codes)
    ci = ci+1;
    
    [sel_trials] = find(VSDI.condition(:,1)==condi);
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    

    back = VSDI.backgr(:,:,VSDI.nonanidx(1));
    
    %to plot single trial
            movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4);
            movie2plot(:,:,end) = back; %clean non-blured background
            
    % -------------------------------------------------------
    % GET MAX-VALUES FOR BOTH METHODS
    % -------------------------------------------------------
    % for 'moviebased_thresh_max'
    tempmax = movmean(movie2plot(:,:,1:end-1),5); %temporal smooth to get the max value par: 3
    
    maxall(ci) = max(tempmax(:));
    
    % For 'wavebased_thresh_local' and 'wavebased_thresh_max'
            dm4_wave = roi_TSave(movie2plot,roimask);
            dm4_wave = movmean(dm4_wave ,5);

    
    wavelocal_max(ci) = max(dm4_wave); 
    
    clear tempmax movie2plot  sel_trials
    
end

maxval = max(maxall);


%----------------------------------------------------------------
... LOOP TO GET TILES, EARLY-PEAK FRAMES AND WAVES
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
    
    back = VSDI.backgr(:,:,VSDI.nonanidx(1));
    
    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4);
    movie2plot(:,:,end) = back; %clean non-blured background
    
    % -------------------------------------------------------
    % CALCULATE dF WAVE FOR dm4 TO ESTABLISH START-END TIMES
    % -------------------------------------------------------
    
    dm4_wave = roi_TSave(movie2plot,roimask);
    dm4_wave = movmean(dm4_wave ,5);
    wave_timebase = VSDI.timebase;
    
    %----------------------------------------------------------------
    % TILES
    %----------------------------------------------------------------
        
    tileset.nrowcol = [1 6];
    tiles.backgr = back;
    tileset.interp =6;
    
    
    fr_pre = 3; %nº of frames before rising the threshold (to set the initial time)
    fr_post = 2;
    
    if plottiles
        % GET REPRESENTATION THRESHOLD ACCORDING TO THE MODE
        %----------------------------------------------------------------
        switch thresh_mode
            case 'moviebased_thresh_max'
                tileset.clims = [0 maxval*fact_clim];
                tileset.thresh = [-maxval*fact_thresh maxval*fact_thresh];
            case 'wavebased_thresh_local'
                localmax = wavelocal_max(ci);
                maxwave = max(wavelocal_max);
                tileset.clims = [0 maxwave*fact_clim]; % the colormap has to be the same  for all conditions to be comparable
                tileset.thresh = [-localmax*fact_thresh localmax*fact_thresh];
            case 'wavebased_thresh_max'
                maxwave = max(wavelocal_max);
                tileset.clims = [0 maxwave*fact_clim]; % the colormap has to be the same for all conditions to be comparable
                tileset.thresh = [-maxwave*fact_thresh maxwave*fact_thresh];
            case 'manual'
                tileset.clims = manual.clims;
                tileset.thresh = manual.thresh;
        end
        
            idx.start = find(dm4_wave > tileset.thresh(2), 1, 'first')- fr_pre;
            idx.end= find(dm4_wave > tileset.thresh(2), 1, 'last') + fr_post;
            
    end    
        
        switch timerange_mode
        case 'refwave'

        tileset.start_ms = VSDI.timebase(idx.start);
        tileset.end_ms = VSDI.timebase(idx.end);
%         tileset.peak_ms = VSDI.timebase(idx.peak);
            case 'manual'
          tileset.start_ms = manual.start_ms; %%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&//////////////////////////77
          tileset.end_ms = manual.end_ms;
        end

        
        if isempty (tileset.start_ms)
%             error('The threshold is too high (lower the value of "fact_thresh"')
            tileset.start_ms = 0;
            tileset.end_ms = 12;
        end
        
        plot_tilemovie_custom_interp2(movie2plot, VSDI.timebase, tileset, [], cmap_tile);
        
        nframes = tileset.nrowcol(1)*tileset.nrowcol(2); %for the title
        
        titulo = [num2str(VSDI.ref) 'cond' num2str(condi) '.clim=' num2str(round(tileset.clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(tileset.thresh(2),1)) '(' num2str(fact_thresh*100) '%)'];
        sgtitle(titulo)
        
        savename= ['tiles' num2str(VSDI.ref) '_' thresh_mode  'trials_factors' num2str(fact_thresh)  '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) '-' num2str(fact_clim) 'pre' num2str(fr_pre) 'post' num2str(fr_post) 'fr' num2str(nframes) movie_ref '.jpg'];
        
        if savetiles
            %                         saveas(gcf, fullfile(savein, savename), 'jpg')
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
            print(fullfile(savein, savename),'-r600','-djpeg') % prints it as you see them
            close
        end
        
        clear idx tileset.start_ms tileset.end_ms
        
    %----------------------------------------------------------------
    % WAVES
    %----------------------------------------------------------------
    if plotwaves
        
        % -------------------------------------------------------
        % CALCULATE %F WAVE FOR EACH ROI
        % -------------------------------------------------------
        
        % DEPRECATED 22/03/22 
        meanF0 = squeeze(mean(VSDmov.F0(:,:,sel_trials),3));
        
        for roii = makeRow(selroi)
            roimask = masks(:,:,roii);
            %                 allroi_waves(:,roii,ci) = roi_TSave(movieave,roimask);
            allroi_waves(:,roii) = roi_TSave_percF_roiwise(movie2plot,roimask, meanF0); %CHANGE TO NOT %F
        end %for roi_i
        
        
        
        % -------------------------------------------------------
        % RANGE
        % -------------------------------------------------------
        
        wave.start_ms = 0; % time in ms for first tile
        wave.end_ms = 1200;
        %           tileset.clims = [-0.9 0.9];
        wave.start= dsearchn(VSDI.timebase, wave.start_ms);
        wave.end = dsearchn(VSDI.timebase, wave.end_ms);
        
        
        % PLOT
        figure
        cmap = roicolors();
        cmap = cmap(1:2:end,:); %roicolors map has double values for 2 hemispheres
        
        back = VSDI.backgr(:,:,VSDI.nonanidx(1));
        
        sp1 = subplot(1,2,1);
        switch roikind
            case 'circle'
                centers = VSDI.roi.circle.center(selroi, :) ;
                roicirc_preview_multiple_cmap(back, centers, VSDI.roi.circle.R, sp1, cmap);
                
            case 'anat'
                roi_preview_multiple(back, VSDI.roi.manual_poly(selroi,:), sp1);
        end
        
        sp1.Visible = 0;
        
        sp2= subplot(1,2,2);
        
        pos1 = get(sp1,'Position');
        pos2 = get(sp2,'Position');
        pos3= [pos2(1) pos2(2) pos1(3) pos1(4)];
        set(sp2, 'Position',pos3)
        
        nroi = length(selroi);
        roicolors= roi_colors();
        
        hold on
        i = 0
        for roii = selroi
            i = i+1;
            waveroi = movmean(allroi_waves(wave.start:wave.end,roii),5);
            plot(VSDI.timebase(wave.start:wave.end), waveroi , 'linewidth', 2, 'Color', cmap(i,:));
            clear waveroi
            %     legend(selroinames{:}, 'Location', 'northeast')
        end
        xlim([wave.start_ms wave.end_ms])
        ylabel('%F')
        
        sgtitle([num2str(VSDI.ref), movie_ref, 'rej', num2str(reject_on), '-cond:' num2str(condi)])
        
        savename= ['waves' num2str(VSDI.ref) movie_ref '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' ];
        
        if savewaves
            saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
            print(fullfile(savein, [savename '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them %STILL TO TEST!
            
            close all
        end
        
        
    end % if plotwaves
    
        
    end
   
    
    
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
% 20/09/22: substantially change the code, delete early, add 'manual' option for timerange 
% 19/09/22: code adjusted so it works also when the threshold is not
% reached (tileset.start_ms is set to zero)
% 17/09/22: add 'manual' option for thresh_mode 
% 12/03/22 : compute_rejectidx function and improve 'cmap'