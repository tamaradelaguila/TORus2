% % ----------------------------------------------------------------
% % % STEP 1:  GET LINES (once for fish and line)
% % %----------------------------------------------------------------
% clearvars -except nfish
% W = pwd;
% cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
% user_settings
% cd(W)
% 
% pathlines = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/raster_fromline';
% 
% load(pathlines)
% 
% nfish = 11;
%     VSDI = TORus('load',nfish);
%     im = VSDI.backgr(:,:,VSDI.nonanidx(1)); axis image
%    imagesc(im); colormap bone
%     roiline = drawline();
%     
%     
%     % GET COORD
%     linecoord(nfish).x = [roiline.Position(1,1) roiline.Position(2,1)]; %x values
%     linecoord(nfish).y = [roiline.Position(1,2) roiline.Position(2,2)]; %y values
%     
%     saveas(gcf, fullfile(pathlines,[num2str(VSDI.ref), 'line.jpg']), 'jpg')
%     close
% 
% save(pathlines, 'linecoord')

% ·······························································
%% STEP 2: GET RASTER FROM (previously drawn) LINE 
% ·······························································
clearvars -except nfish
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
pathlines = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/raster_fromline';

load(pathlines)

saveraster =0;
savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/raster_fromline';


nfish = 11;
% fact_thresh =0.4; % @SET : limits parameters
fact_clim= 1.2;

cond_codes = [400:404];

movie_ref = '_18filt6'; %

%----------------------------------------------------------------
% @SET: REJECT SETTINGS 
%----------------------------------------------------------------
reject_on= 3;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1 ; %
setting.GSabsthres_reject = 1; %
setting.forcein = 1; %
setting.bradyvisual = 0;

%% LOAD / COMPUTE SETTINGS
%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus('load', nfish);
VSDmov = TORus('loadmovie',nfish,movie_ref);

%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------
rej = 'reject' ;
if reject_on > 1
rej = [rej num2str(reject_on)];
end

rejectidx = [];
reject_idx  = compute_rejectidx(VSDI, reject_on, setting);

%----------------------------------------------------------------------
% GET MAX FROM ALL CONDITIONS TO SET AS REPRESENTATION'S THRESHOLD AND
% COLORLIMITS
%----------------------------------------------------------------------
ci = 0;
for condi = makeRow(cond_codes)
    ci = ci+1;
    
    [sel_trials] = find(VSDI.condition(:,1)==condi);
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    
    
    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4);
    back = VSDI.backgr(:,:,VSDI.nonanidx(1));
    movie2plot(:,:,end) = back; %clean non-blured background after calculating the F0 value (jsut for visualization purposes)

    % SMOOTH DATA TO GET THE MAX VALUE
    tempmax = movmean(movie2plot(:,:,1:end-1),5); %temporal smooth to get the max value par: 3
    
    maxall(ci) =  max(tempmax(:));
    
end

maxval = max(maxall);
trange = 1:270;
    %% ----------------------------------------------------------------
    % GET  RASTER
    %----------------------------------------------------------------

    for condi = makeRow(cond_codes)

    % APPLY TRIALS REJECTION 
    %----------------------------------------------------------------
    [sel_trials] = find(VSDI.condition(:,1)==condi);
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    
    % GET AVERAGE MOVIE 
    %----------------------------------------------------------------
    movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4);
    back = mean(VSDmov.data(:,:,end,sel_trials),4);
    

    % GET COORD FROM SAVED STRUCTURE
    %----------------------------------------------------------------
    xs = linecoord(nfish).x;
    ys = linecoord(nfish).y;
    
%     F0val = improfile(back, xs, ys);

    % GET (%F0) PIXEL VALUES FOR EACH TIMEPOINT AND STORE IN RASTER
    %----------------------------------------------------------------
    for ti = 1:length(VSDI.timebase) %leave background out
        im = movie2plot(:,:,ti);
        pixval = improfile(im, xs, ys);
%         rasterline(:,ti) = pixval./F0val;
        rasterline(:,ti) = pixval;
       clear im
    end
    
    % PLOT
    %----------------------------------------------------------------
    clim = [0 maxval*fact_clim];
%     clim = [0 14];

    figure
    % plot(pixval, 1:length(pixval)); axis tight; title('pixel values')
    ax1= subplot(1,2,1);
    imagesc(VSDI.timebase(trange), [], rasterline(:,trange));
    xline(0, '--w', 'Linewidth', 1.3)
    colormap(ax1, parula)
    colorbar;
    set(ax1, 'clim', clim)
    ylabel('caudo-rostral pixels')
    xlabel('s')
    title('raster of activity')

    ax2= subplot(1,2,2);
    
    coord = [xs(1) ys(1); xs(2) ys(2)];
    im = VSDI.backgr(:,:,VSDI.nonanidx(1));
%     im = imrotate(im,90)
    imagesc(im);
    axis image
    ax2 = gca ;
    colormap(ax2, 'bone')
    drawline(ax2, 'Position', coord, 'Color', [1 1 1])
    title('pixels selected')
        
        
        sgtitle([num2str(VSDI.ref), movie_ref, '.rej', num2str(reject_on), '.c = ' num2str(condi) '(' num2str(VSDI.condition(sel_trials(1),4)) 'mA)'])
        
        savename= ['rasterline' num2str(VSDI.ref) movie_ref num2str(condi) '_rej'  num2str(reject_on)];
        
        if saveraster
            saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
%             set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%             print(fullfile(savein, [savename '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them %STILL TO TEST!
            close all
        end
        
end %condi

blob()
blob() ; pause(0.2) ; blob()

% newsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes_code/03_figure_sketch/def_figs/tiles';
% saveas(gcf, fullfile(newsave,savename) 'jpg')

%-----------------------------------------------------------------
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


% Created: 10/02/22
% source: /home/tamara/Documents/MATLAB/VSDI/seltrials/plot_code/informes_code/03_figure_sketch/plot_tiles_and_waves.m
% Last update: 