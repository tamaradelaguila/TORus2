%% PLOT MOVIES WITH
clear
user_settings
%
% roi2draw = VSDI.roi.manual_poly([1,2,7,8,13,14]);
% plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset, roi2draw);
% plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
%
% ref_labels0 = {'line1', 'line2'}
%  [coord0] = draw_brainref(VSDI.crop.preview,ref_labels0)

%% OVERLAY ANY NUMBER OF TILES : from a movie
clearvars -except path
path2save = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/tiles';

nfish = 1;
VSDI = TORus ('load', nfish)
% movie_ref = '_12filt5'; % input movie
% movie_ref = '_13filt6'; %
% movie_ref = '_14filt6';
% movie_ref = '_15filt5'; %

movie_ref = '_17filt5'; %@ SET

VSDmov = TORus('loadmovie',nfish,movie_ref);

%@ SET
% selectmode = 'all';
% selectmode = 'nosharks';
selectmode = 'onlysharks';

saveandclose = 0; %@ SET
%% SELECT CASES

trialkind = unique(VSDI.condition(:,1));
trialkind = makeRow(trialkind(~isnan(trialkind)));


%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

reject_on = [0];  %@ SET
% Subsettings:
setting.manual_reject = 0; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET+
setting.force_include = 0; %@ SET

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

for condi = trialkind %to make sure that only the conditions of interest are computed)
    condtrials = makeCol(find(VSDI.condition(:,1)==condi));
    sel_trials  = [sel_trials; condtrials];
end
sel_trials= sort(sel_trials);

switch selectmode
    case 'all'
        sel_trials = sel_trials;
        
    case 'nosharks'
        sel_trials = makeRow(setdiff(sel_trials, rejectidx));
        
    case 'onlysharks'
        sel_trials = makeRow(VSDI.reject.GSabs025);
        
end

%% DRAW ROIS BEFORE PLOTTING
% refer = {'1', '2'}'; %set a name for each reference
%
% [roi2draw, ~] = roi_draw(VSDI.crop.preview,refer); %in each column, one set of rois

% %% PLOT AND SAVE LOOP
%
% % -----------------------------------------
% % 'plot_tilemovie_custom' SETTINGS
% % -----------------------------------------
%
%     tileset.start_ms = 0; % time in ms for first tile
%     tileset.end_ms = 250;
%     %           tileset.clims = [-0.9 0.9];
%     tileset.clims = [-6 6];
%     tileset.thresh = [-1.7 1.7];
%
%     tileset.nrowcol = [4 4];
%     ntiles = tileset.nrowcol(1)*tileset.nrowcol(2);
%     tileset.backgr = VSDI.crop.preview; %use non-meaned background
%
% % PLOT AND SAVE TILES
%
% for condi =   makeRow(trialkind) % [1000:1003]
%     [idxA] = find(VSDI.condition(:,1)==condi);
%     idxA = intersect(idxA, sel_trials);
%
%     % cond0 = ();
%     % [idx0] = find(VSDI.condition(:,1)==400);
%
%     %to plot single trial
%     movie2plot = mean(VSDmov.data(:,:,:,idxA),4) ;
%
%     movie2plot_dif = mean(VSDmov.data(:,:,:,idxA),4); % - mean(VSDmov.data(:,:,:,idx0),4) ;
%
%     % PLOT
% %     roi2draw = VSDI.roi.manual_poly([1,2,7,8,13,14]);
% %     plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset, roi2draw);
%     plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
%
% %         plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
%
% %             plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
%
%
% %     titulo = [num2str(VSDI.ref) '(' movie_ref(2:end)  ') cond' num2str(condi) ' clim=' num2str(tileset.clims(1)) ';thresh=' num2str(tileset.thresh(1)) selectmode '-trials n=' num2str(length(idxA))];
% %     sgtitle(titulo)
%
%     name2save = ['TILES_' num2str(VSDI.ref) movie_ref '_' selectmode  '_cond' num2str(condi) '(n' num2str(length(idxA)) ')_clim' num2str(tileset.clims(1)) 'thresh' num2str(tileset.thresh(1)) '_ntiles' num2str(ntiles)];
%     saveas(gcf,fullfile(path2save ,name2save),'jpg')
%     close
%
% end
%
% blob()

% %% PLOT 1 ROW + WAVES
%
% % PLOT-FUNCTION SETTINGS
% roi1 = 1;
% roi2 = 7;
%
% id1= roi1;
% id2= roi2;
% movieset.times = VSDI.timebase;
% movieset.backgr = VSDI.crop.preview;
%
% tileset.start_ms = 0;
% tileset.end_ms = 250;
% tileset.time2plot = 0; %select time (ms)
% %         tileset.clims = [-0.2 0.2];%  values above or below will be saturated
% %         tileset.thresh = [-0.08 0.08]; %% values inside this range will be zero-out
% tileset.ntiles = 4;
% % %
% tileset.clims = [-12 12];% LOWER values above or below will be saturated
% tileset.thresh = [-4 4]; %%LOWER values inside this range will be zero-out
% % %
% waveset.coord =  VSDI.roi.circle.center([id1 id2],:);
% waveset.r = VSDI.roi.circle.R;
% waveset.xlim = [0 600];
% waveset.ylim =tileset.clims;
%
% for condi =  makeRow(trialkind) % [1000:1003]
%     [idxA] = find(VSDI.condition(:,1)==condi);
%     idxA = intersect(idxA, sel_trials);
%
%     % cond0 = ();
%     % [idx0] = find(VSDI.condition(:,1)==400);
%
%     %to plot single trial
%     movie2plot = mean(VSDmov.data(:,:,:,idxA),4) ;
%
%
%     movieset.data =movie2plot;
%
%     % END OF FUNCTION SETTINGS
%
%     % ------------------
%     % APPLY FUNCTION
%     % ------------------
%     tilesrow_nframes_nwaves(movieset, tileset,waveset)
%     % ------------------
%
%     % ADJUSTMENTS FOR THE PLOT
%     % Plot the threshold and clims:
%     yline(tileset.thresh(2), 'b--')
%     yline(tileset.clims(2), 'r--')
%
%     tempmA = num2str(VSDI.condition(idxA(1),4));
%     titulo = [num2str(VSDI.ref), ':', tempmA,'mA (cod=',num2str(condi), ')' ];
%     sgtitle(titulo)
%     set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 3]);
%
%
%     name2save = ['TILES_waves' num2str(VSDI.ref) movie_ref '_' selectmode  '_cond' num2str(condi) '(n' num2str(length(idxA)) ')_clim' num2str(tileset.clims(1)) 'thresh' num2str(tileset.thresh(1)) '_ntiles' num2str(tileset.ntiles)];
%     saveas(gcf,fullfile(path2save ,name2save),'jpg')
%     close
%
% end
%
% blob()

%% BOTH TILES GRID + {TILES ROW & WAVES}

% -----------------------------------------
% TILES GRID ('plot_tilemovie_custom') SETTINGS
% -----------------------------------------

tileset.start_ms = 0; % time in ms for first tile
tileset.end_ms = 700;
%           tileset.clims = [-0.9 0.9];

climabs =  8; %@ SET
thresabs = 1.5; %@ SET

tileset.clims = [-climabs climabs];
tileset.thresh = [-thresabs thresabs];

tileset.nrowcol = [6 5];
ntiles = tileset.nrowcol(1)*tileset.nrowcol(2);
tileset.backgr = VSDI.crop.preview; %use non-meaned background


% -----------------------------------------
% TILES ROW ('tilesrow_nframes_nwaves') SETTINGS
% -----------------------------------------
roi1 = 3;
roi2 = 7;

id1= roi1;
id2= roi2;
movieset.times = VSDI.timebase;
movieset.backgr = VSDI.crop.preview;

tileset2.start_ms = 0;
tileset2.end_ms = 130;
tileset2.time2plot = 0; %select time (ms)
%         tileset2.clims = [-0.2 0.2];%  values above or below will be saturated
%         tileset2.thresh = [-0.08 0.08]; %% values inside this range will be zero-out
tileset2.ntiles = 4;
% %
tileset2.clims = tileset.clims
tileset2.thresh = tileset.thresh
% %
waveset.coord =  VSDI.roi.circle.center([id1 id2],:);
waveset.r = VSDI.roi.circle.R;
waveset.xlim = [0 600];
waveset.ylim = tileset.clims;


% PLOT AND SAVE TILES

for condi =   404 %makeRow(trialkind) % [1000:1003]
    [idxA] = find(VSDI.condition(:,1)==condi);
    idxA = intersect(idxA, sel_trials);
    
    % cond0 = ();
    % [idx0] = find(VSDI.condition(:,1)==400);
    
    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,idxA),4) ;
    
    
    movieset.data =movie2plot;
    
    % END OF FUNCTION SETTINGS
    
    % ------------------------
    % TILES GRID
    % -------------------------
    tilesrow_nframes_nwaves(movieset, tileset2,waveset)
    % ------------------
    
    % ADJUSTMENTS FOR THE PLOT
    % Plot the threshold and clims:
    yline(tileset.thresh(2), 'b--')
    yline(tileset.clims(2), 'r--')
    
    tempmA = num2str(VSDI.condition(idxA(1),4));
    titulo = [num2str(VSDI.ref), ':', tempmA,'mA (cod=',num2str(condi), ')' ];
    sgtitle(titulo)
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 3]);
    
    
    if saveandclose
    name2save = ['TILES_waves' num2str(VSDI.ref) movie_ref '_' selectmode  '_clim' num2str(tileset2.clims(1)) 'thresh' num2str(tileset2.thresh(1)) '_ntiles' num2str(tileset2.ntiles)  '_cond' num2str(condi) '(n' num2str(length(idxA)) ')' '.jpg'];
    saveas(gcf,fullfile(path2save ,name2save),'jpg')
    close
    end
    
    % ------------------------
    % TILES GRID
    % -------------------------
    %     roi2draw = VSDI.roi.manual_poly([1,2,7,8,13,14]);
    %     plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset, roi2draw);
    plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
    
    %         plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
    
    %             plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
    
    
    %     titulo = [num2str(VSDI.ref) '(' movie_ref(2:end)  ') cond' num2str(condi) ' clim=' num2str(tileset.clims(1)) ';thresh=' num2str(tileset.thresh(1)) selectmode '-trials n=' num2str(length(idxA))];
    %     sgtitle(titulo)
    
        if saveandclose
    name2save = ['TILES_' num2str(VSDI.ref) movie_ref '_' selectmode  '_clim' num2str(tileset.clims(1)) 'thresh' num2str(tileset.thresh(1)) '_ntiles' num2str(ntiles)  '_cond' num2str(condi) '(n' num2str(length(idxA)) ')' '.jpg'];
    disp(name2save)
    saveas(gcf,fullfile(path2save ,name2save),'jpg')
    close
        end
    
end

blob()

% updated: 26/10/21