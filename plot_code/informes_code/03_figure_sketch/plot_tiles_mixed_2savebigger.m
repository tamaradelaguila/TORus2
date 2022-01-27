clear
user_settings

%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
nfish = 8;
cond_codes = [1003];

plottiles = 0;
savetiles = 0;
plotwaves=1;
savewaves = 0;

movie_ref = '_18filt6'; % input movie '_17filt5'
% movie_ref = '_17filt5'; % input movie '_17filt5'

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/def_figs/tiles';

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on= 1;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.force_include = 0; %

%----------------------------------------------------------------
% @SET: for roi-waves plot
%----------------------------------------------------------------
% selroinames = {'dm4m_R',  'dm2_R', 'dm3_R' ,'dm1_R','dldm_R'};%dm3,
selroinames = {'dm4m_R2',  'dm2_R2' ,'dm1_R','dldm_R2'};%dm3,
roikind = 'circle'; %
% roikind = 'anat';

%% LOAD / COMPUTE SETTINGS
%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus ('load', nfish);
VSDmov = TORus('loadmovie',nfish,movie_ref);

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
end

%----------------------------------------------------------------
%% CODE: TILES + WAVEs (all roi for each condition)
%----------------------------------------------------------------
% cond_codes = [201];
fact_thresh =0.4; % @SET : limits parameters
fact_clim= 1.2;
        

for condi = cond_codes
    
    %----------------------------------------------------------------
    % 
    %----------------------------------------------------------------
    
    [sel_trials] = find(VSDI.condition(:,1)==condi);
    sel_trials = setdiff(sel_trials, rejectidx);
    
    back = VSDI.backgr(:,:,VSDI.nonanidx(1));

    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4);
    movie2plot(:,:,end) = back; %clean non-blured background
    
    
    %----------------------------------------------------------------
    % TILES
    %----------------------------------------------------------------
    % settings
    tileset.start_ms = 0; % time in ms for first tile
    tileset.end_ms = 700;
    %           tileset.clims = [-0.9 0.9];
    idx.start= dsearchn(VSDI.timebase, tileset.start_ms);
    idx.end = dsearchn(VSDI.timebase, tileset.end_ms);

    
    if plottiles
        
        
        
        % SMOOTH DATA TO GET THE MAX VALUE
        maxval = movmedian(movie2plot(:,:,1:end-1),8 ,3); %temporal smooth to get the max value
        maxval = max(maxval(:));
        
        
        tileset.clims = [-maxval*fact_clim maxval*fact_clim];
        %     tileset.clims = [-8 8];
        %     tileset.thresh = [-maxval*0.4 maxval*0.4];
        tileset.thresh = [-maxval*fact_thresh maxval*fact_thresh];
        %     tileset.thresh = [-4.5 4.5];
        
        tileset.nrowcol = [1 6];
        tiles.backgr = back;
        tileset.interp =6;
        
        % get info to plot in title
        
        % plot
        %     roi2draw = VSDI.roi.manual_poly([1,2,7,8,13,14]);
        %     plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset, roi2draw);
        plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
        
        %         plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
        
        %             plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
        
        
        titulo = [num2str(VSDI.ref) 'cond' num2str(condi) '.clim=' num2str(round(tileset.clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(tileset.thresh(2),1)) '(' num2str(fact_thresh) '%)'];
        sgtitle(titulo)
        
        savename= ['tiles' num2str(VSDI.ref) movie_ref 'cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_factors' num2str(fact_thresh) '-' num2str(fact_clim) '.jpg'];
        
        if savetiles
            saveas(gcf, fullfile(savein, savename), 'jpg')
            close
        end
        
    end
    
    %----------------------------------------------------------------
    % WAVES
    %----------------------------------------------------------------
    if plotwaves
        
            % -------------------------------------------------------
            % CALCULATE WAVE FOR EACH ROI 
            % -------------------------------------------------------
            meanF0 = squeeze(mean(VSDmov.F0(:,:,sel_trials),3));

            for roii = makeRow(selroi)
                roimask = masks(:,:,roii);
%                 allroi_waves(:,roii,ci) = roi_TSave(movieave,roimask);
                allroi_waves(:,roii) = roi_TSave_percF_roiwise(movie2plot,roimask, meanF0);
            end %for roi_i
            
            % -------------------------------------------------------
            % CALCULATE WAVE FOR EACH ROI
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
    i = 1;
    for roii = selroi
    plot(VSDI.timebase(wave.start:wave.end),allroi_waves(wave.start:wave.end,roii), 'linewidth', 1.3, 'Color', cmap(i,:));
    i = i+1; %in this way the colours of the roi are respected
%     legend(selroinames{:}, 'Location', 'northeast')
    end
    
    ylabel('%F')
 
            sgtitle([num2str(VSDI.ref), movie_ref, 'rej', num2str(reject_on), '-cond:' num2str(condi)])
            
            savename= ['waves' num2str(VSDI.ref) movie_ref 'cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials.jpg'];

            if savewaves
            saveas(gcf, fullfile(savein, savename), 'jpg')
            close all
            end


    end % if plotwaves
    
end % for condi


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

