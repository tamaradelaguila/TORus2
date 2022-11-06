close all
clearvars -except nfish
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)
%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------

nfish = 12;
cond_codes = [400:404];

fact_clim= 1;

% times to plot
t.ini = -300; %ms
t.end = 800;

movie_ref = '_18filt6';
saveraster =1;
name_ending = 'dm4_dm2_R_smaller'; % ATT that matches to the intention

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/raster_fromline';

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------
reject_on= 3;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1 ; %
setting.GSabsthres_reject = 1; %
setting.forcein = 1; %
setting.bradyvisual = 0;

%----------------------------------------------------------------
%% STEP 0: LOAD DATA
%----------------------------------------------------------------
VSDI = TORus('load', nfish);
VSDmov = TORus('loadmovie',nfish,movie_ref);

%% ----------------------------------------------------------------
% % STEP 1:  GET LINES (once for fish and line)
% %----------------------------------------------------------------
im = VSDI.backgr(:,:,VSDI.nonanidx(1)); axis image
imagesc(im); colormap bone

for n = 1:2
title(['draw line nº' num2str(n) '(' name_ending ')'], 'Interpreter', 'none')
roiline = drawline();
pause

% GET COORD
linecoord(n).x = [round(roiline.Position(1,1)) round(roiline.Position(2,1))]; %x values
linecoord(n).y = [round(roiline.Position(1,2)) round(roiline.Position(2,2))]; %y values
end

% ·······························································
%% STEP 2: GET RASTER FROM (previously drawn) LINE
% ·······························································


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

idx.ini = dsearchn(VSDI.timebase, t.ini) ;
idx.end= dsearchn(VSDI.timebase, t.end) ;
trange = idx.ini:idx.end;
%% ----------------------------------------------------------------
% GET  RASTER
%----------------------------------------------------------------


for condi = makeRow(cond_codes)
    
        % GET SAVED COORDINATES
        %----------------------------------------------------------------
        
        x1 = linecoord(1).x;
        y1 = linecoord(1).y;

        x2 = linecoord(2).x;
        y2 = linecoord(2).y;
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
    %     F0val = improfile(back, xs, ys);
    
    % GET (%F0) PIXEL VALUES FOR EACH TIMEPOINT AND STORE IN RASTER
    %----------------------------------------------------------------
    for ti = 1:length(VSDI.timebase) %leave background out
        im = movie2plot(:,:,ti);
        %         rasterline(:,ti) = pixval./F0val;
        rasterline1(:,ti,condi) = improfile(im, x1, y1);
        rasterline2(:,ti,condi) = improfile(im, x2, y2);

        clear im pixval
    end
    
end %for condi

%% ----------------------------------------------------------------
% PLOT
%------------------------------------------------------------------
close all

temp1 = rasterline1(:,trange,:);
maxval1 = max(temp1(:));
temp2 = rasterline2(:,trange,:);
maxval2 = max(temp2(:));

maxval = max([maxval1, maxval2]);

clim = [0 maxval*fact_clim];
cmap = lines;

%     clim = [0 14];
for condi = makeRow(cond_codes)
    figure
    % plot(pixval, 1:length(pixval)); axis tight; title('pixel values')
    ax1= subplot(3,1,1);
    imagesc(VSDI.timebase(trange), [], rasterline1(:,trange, condi));
    xline(0, '--w', 'Linewidth', 1.3)
    colormap(ax1, parula)
    colorbar;
    set(ax1, 'clim', clim)
    ylabel('pixels', 'Color', cmap(3,:))
    xlabel('s')
    title('raster of activity','Color', cmap(3,:))
    
    ax2= subplot(3,1,2);
    imagesc(VSDI.timebase(trange), [], rasterline2(:,trange, condi));
    xline(0, '--w', 'Linewidth', 1.3)
    colormap(ax2, parula)
    colorbar;
    set(ax2, 'clim', clim)
    ylabel('pixels', 'Color', cmap(4,:))
    xlabel('s')
    title('raster of activity','Color', cmap(4,:))

    
    ax3= subplot(3,1,3);
    
    im = VSDI.backgr(:,:,VSDI.nonanidx(1));
    %     im = imrotate(im,90)
    imagesc(im);
    axis image
    ax3 = gca ;
    colormap(ax3, 'bone')
    for n = 1:2
        xs = linecoord(n).x;
        ys = linecoord(n).y;
%         coord = [xs(1) ys(1); xs(2) ys(2)];
%         drawline(ax3, 'Position', coord1, 'Color',cmap(n+2,:))
     line(xs,ys, 'Linewidth', 2 ,'Color', cmap(n+2,:) )
    end

    title('pixels selected')
    
    sgtitle([num2str(VSDI.ref), movie_ref, '.rej', num2str(reject_on), '.c = ' num2str(condi) '(' num2str(VSDI.condition(sel_trials(1),4)) 'mA)'])
    
    savename= ['rasterline' num2str(VSDI.ref) movie_ref num2str(condi) '_rej'  num2str(reject_on) '_x=' num2str(linecoord.x) 'y=' num2str(linecoord.y) 'Fclim' num2str(fact_clim)];
    
    if saveraster
        saveas(gcf, fullfile(savein, [savename '_' name_ending '.jpg']), 'jpg')
        %             set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
        %             print(fullfile(savein, [savename '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them %STILL TO TEST!
        close all
    end
    
end %for condi

blob()

return
%% PRINT 2 REFERENCE LINES FROM MANUAL COORD
x1 = [36 40];
y1 = [8 17];

x2 = [30 49];
y2 = [12 23];

cmap = lines;
%     coord1 = [x1(1) y1(1); x1(2) y1(2)];
%     coord2 = [x2(1) y2(1); x2(2) y2(2)];
im = VSDI.backgr(:,:,VSDI.nonanidx(1));
%     im = imrotate(im,90)
imagesc(im);
axis image
ax = gca ;
colormap(ax, 'bone')
%     drawline(ax, 'Position', coord1, 'Color', cmap(3,:))
%         drawline(ax, 'Position', coord2, 'Color', cmap(4,:))
line(x1,y1, 'Linewidth', 3 ,'Color', cmap(3,:) )
line(x2,y2, 'Linewidth', 3 ,'Color', cmap(4,:) )

title('pixels selected')


%% Created: 11/02/22
% source: /home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/03_figure_sketch/raster_from_savedline.m
% Last update: