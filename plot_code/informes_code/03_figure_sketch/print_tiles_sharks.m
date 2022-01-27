%% PRINT SHARKS TILES (AVERAGE MOVIE)
clear

%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
nfish = 1;
cond_codes = [302:303];
% message = 'tren largo - no aleat (intra?)';
message = 'sharks';

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on= 0;

setting.manual_reject = 0; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.force_include = 0; %

%----------------------------------------------------------------
% LOAD PARAMETERS
%----------------------------------------------------------------
VSDI = TORus ('load', nfish);
movie_ref = '_17filt5'; % input movie
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
% PLOT TILES 
%----------------------------------------------------------------
for condi = cond_codes
        
%%    
    [sel_trials] = find(VSDI.condition(:,1)==condi);
    
%     local_reject = intersect(sel_trials, rejectidx);%just to display
%     later; MUTED
%     sel_trials = setdiff(sel_trials, rejectidx); MUTED

    seltrials = intersect(sel_trials, makeRow(VSDI.reject.GSabs025)); %%% ADAPTED

    
    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4);
    movie2plot(:,:,end) = VSDI.crop.preview; %clean non-blured background
    % settings
    tileset.start_ms = 726; % time in ms for first tile
    tileset.end_ms = 1434;
    tileset.clims = [-25 25];
    tileset.thresh = [-2 2];
    
    tileset.nrowcol = [10 10]; % [31 10];
    tiles.backgr = VSDI.crop.preview;
    % get info to plot in title
    
    % plot
%     roi2draw = VSDI.roi.manual_poly([1,2,7,8,13,14]);
%     plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset, roi2draw);
    plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);
%     tileset.timeidx = 
%     devo_plot_allframes(movie2plot, tileset);

%         plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);

%             plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);

        
    titulo = [num2str(VSDI.ref) 'cond' num2str(condi)];
    sgtitle(titulo)
   
    name = [titulo '_' num2str(tileset.start_ms) ' to' num2str(tileset.end_ms) '.jpg'];
    saveas(gcf, name, 'jpg')
    close 
end


% plot_allframes(movie2plot,tileset)
%----------------------------------------------------------------
% PRINT INCLUDED AND EXCLUDED TRIALS
% ----------------------------------------------------------------

for condi = cond_codes
    disp(['Included trials for condition' num2str(condi) ':'])
    disp(sel_trials)
    disp('%')
    disp(['Rejected trials for condition' num2str(condi) ':'])
    disp(local_reject)
end

%% ----------------------------------------------------------------
%----------------------------------------------------------------
% DISPLAY PARAMETERS
% ----------------------------------------------------------------
%----------------------------------------------------------------
% clearvars -except measure reject_on nfish cond_codes VSDI movies waves ref_wave ...
%     window window_ave noise noise_ave rejectidx roiname1 roiname2 method lat_limit...
%     wavesylim selroi setting
display(['reject' num2str(reject_on)])

disp('%%%SETTINGS')
disp(movie_ref)
disp(['clim:' num2str(tileset.clims)])
disp(['thresh:' num2str(tileset.thresh)])
%% 
% namepdf = [num2str(VSDI.ref) '-block' num2str(cond_codes(1)) '.pdf'];
% % publish('source_code_mixed_pdf.mlx','pdf');
%
% mlxloc = strcat(pwd,'/source_code_mixed_pdf.mlx');
% fileout = strcat(pwd,namepdf);
% matlab.internal.liveeditor.openAndConvert(mlxloc,fileout);
blob()
