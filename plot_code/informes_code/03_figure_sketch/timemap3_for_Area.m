clear
close all
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

%----------------------------------------------------------------
% 'TIMESPAN MODE'
%----------------------------------------------------------------

% 'method_1': to make sure you have both regions
...'t0' when the most active reaches the threash;
    ...'tend' when the less active crosses it
    ... makes sense with thresh_mode = 'smallerwave'
    
    % 'method_2'(~tiles): relative to most active's start and end
...'t0' when the most active reaches the threash;
    ...'tend' when the most active crosses the threh again
    ... makes sense with thresh_mode = 'biggerwave'
    
% 'method_3': relative to most active's start and peak
...'t0' when the most active reaches the threash;
    ...'tend' when the most active reaches the peaktime 
    
% 'fixed_timerange' - manually introduce ms for start and end times

% 'THRESH MODE'
% 'maxmovie'
...'thresh set to the max of the (temporally smoothed) movie (as done with the tiles initially)'
                    
% 'biggerwave'
...the largest wave is the reference wave respect to which the thresh is set 
    
%'smallerwave'
...the smallest wave is the reference wave respect to which the thresh is set 

%----------------------------------------------------------------

%----------------------------------------------------------------
% @SET: fish 
%----------------------------------------------------------------
nfish = 11;
dataunits = '%dF'; % '%dF' 'dF'

%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus ('load', nfish);

%% 
close all
clearvars -except path nfish  VSDI VSDmov ref_movie dataunits

%----------------------------------------------------------------
% @SET: basic settings
%----------------------------------------------------------------

condi = 402; % DO NOT INCLUDE THE BLANK CONDITION (or the threshold will fail)

save_timemap = 0;
highdefinition = 1;

%----------------------------------------------------------------
% @SET: CONTOURS PARAMETERS
%----------------------------------------------------------------
roiC_names = {'dm4m_R2', 'dm4m_L2', 'dm2_R2', 'dm2_L2'}; % reference rois for CONTOURS determination.
... the most and less active area will be determined among these
roiA_names = {'dm4_R', 'dm4_L', 'dm2_R', 'dm2_L'}; % for AREA determination (those with which intersect the active pixels)

nframes = 6;

fact_thresh =0.96; % @SET : limits parameters
thresh_mode = 'smallerwave'; % 'maxmovie', 'biggerwave', 'smallerwave'

timespan_mode =  'fixed_timerange' ;  %  'fixed_timerange' 'method_1' 'method_2' 'method_3'
manual_timespan = [108  204]; % (ms) will be used only in 'fixed_timerange' timespan_mode

% DETAILS  IN HEADER

% nº of additional frames before and after startime and endtime,
% respectively
fr_pre = 3 ; %nº of frames before rising the threshold (to set the initial time)
fr_post = 2;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/def_figs/timemaps2_area/'; % CAHNGEEEEEE

% DETERMINE MOVIES REFERENCES ACCORDING TO THE UNITS USED
switch dataunits
    case 'dF'
        ref_movie = '_18filt6'; % input movie '_17filt5'
    case '%dF'
        % ref_movie = '_21filt6';
        % ref_movie = '_22filt7';
        % ref_movie = '_23filt8';
        ref_movie =  '_24filt9';
end

VSDmov = TORus('loadmovie',nfish,ref_movie);

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on= 4;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.forcein = 0; %


% COMPUTE SETTINGS

%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------

reject_idx  = compute_rejectidx(VSDI, reject_on, setting);

%----------------------------------------------------------------
% SET ROI REFERENCES AND SELECT ROI
%----------------------------------------------------------------

roiC_labels = VSDI.roi.labels_circ;
roiA_labels = VSDI.roi.labels;

roiC_mask = VSDI.roi.circle.mask;
roiA_mask = VSDI.roi.manual_mask; 

roiC_idx = name2idx(roiC_names, roiC_labels);
roiA_idx = name2idx(roiA_names, roiA_labels);

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

%----------------------------------------------------------------
... GET AVERAGE MOVIE
    %----------------------------------------------------------------
[sel_trials] = find(VSDI.condition(:,1)==condi);

if reject_on
    sel_trials = setdiff(sel_trials, reject_idx);
end

back = VSDI.backgr(:,:,VSDI.nonanidx(1));

%to plot single trial
movieave = mean(VSDmov.data(:,:,:,sel_trials),4);
movieave(:,:,end) = back; %clean non-blured background


%----------------------------------------------------------------
% GET THRESHOLD
%----------------------------------------------------------------
for ri = 1:numel(roiC_idx)
    roi_idx = roiC_idx(ri);
    roi_mask = roiC_mask(:,:,roi_idx);
    
    % GET WAVE AND SMOOTH DATA
    tempwave = roi_TSave(movieave,roi_mask);
    wave(:,ri) = movmean(tempwave,5); %temporal smooth to get the max value par: 3
    
    tempR = devo_peak2peak(wave(:,ri), VSDI.timebase , feedf.window, [], feedf.method, 0, 0);
    
    maxall(ri) = max(wave(:,ri));
    peakidx_all(roi_idx) = tempR.peakidx(2);
    
end

[minvalue, minidx] = min(maxall); % the values need to be alone in a vector, instead of a vector respecting the roi's real idx because otherwise the minimum would be set to zero
[maxvalue, maxidx] = max(maxall); % 

ref_small = roiC_names{minidx} ; % REZPECT TO WHICH SET THE THRESHOLD
ref_big = roiC_names{maxidx} ; %

% --------------------------------------------------------------------------------
%  SET START-END TIMES
% --------------------------------------------------------------------------------

% roiC_names{4}
idx_roismall = name2idx(ref_small, roiC_labels);
idx_roibig = name2idx(ref_big, roiC_labels);

% GET MASKS
mask_small = roiC_mask(:,:,idx_roismall);
mask_big = roiC_mask(:,:,idx_roibig);

% GET WAVES
wavesmall = roi_TSave(movieave,mask_small);
wavesmall = movmean(wavesmall ,5);

wavebig = roi_TSave(movieave,mask_big);
wavebig = movmean(wavebig ,5);

wave_timebase = VSDI.timebase;

fr_pre = 0; %nº of frames before rising the threshold (to set the initial time)
fr_post = 2;

switch thresh_mode
    case  'maxmovie'
        tempmax = movmean(movieave(:,:,1:end-1),5); %temporal smooth to get the max value par: 3
        moviemax= max(tempmax(:));
        disp('thresh set to the max of the (temporally smoothed) movie (as in tiles)')
            
        thresh = moviemax*fact_thresh ;
        
    case 'biggerwave'
    thresh = maxvalue*fact_thresh ;
    disp(['reference wave to set thresh:' ref_big])
    
    case 'smallerwave'
        thresh = minvalue*fact_thresh ;
        disp(['reference wave to set thresh:' ref_small])
end

switch timespan_mode
    
            
    case   'fixed_timerange'
        
        idx.start = dsearchn(VSDI.timebase, manual_timespan(1));
        idx.end = dsearchn(VSDI.timebase, manual_timespan(end));

    case 'method_1'
        disp([ 'Ref to set start:' ref_big '. ref to set end time:' ref_small])
        
        idx.start = find(wavebig > thresh, 1, 'first') - fr_pre;
        idx.end= find(wavesmall > thresh, 1, 'first') + fr_post;
        
        if isempty (idx.end)
            error(['The threshold is too high for' ref_peak])
        end
        
    case 'method_2'
        disp(['Ref to set start:' ref_big '. ref to set end time:' ref_small])
        
        idx.start = find(wavebig > thresh, 1, 'first') - fr_pre;
        idx.end = find(wavebig > thresh, 1, 'last') + fr_post;
        
    case 'method_3'
        idx.start = find(wavebig > thresh, 1, 'first') - fr_pre;
        idx.end = peakidx_all(idx_roibig);
end

%     titulo = [num2str(VSDI.ref) 'cond' num2str(condi) '.clim=' num2str(round(tileset.clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(thresh(2),1)) '(' num2str(fact_thresh) '%)'];
%     sgtitle(titulo)

%----------------------------------------------------------------
... GET CONTOURS AND INTERSECTION AREA
    %------------------------------------------------------------
frameidx = round( linspace(idx.start,idx.end, nframes));

backgr = movieave(:,:,end);
if highdefinition
    backgr = interp2(backgr,5, 'cubic');
end

ti = 0; % counter for time
for fri = 1:numel(frameidx)
    
    frame = movieave(:,:,frameidx(fri));
    
    if highdefinition
        frame = interp2(frame,5, 'cubic');
    end
    
    
    % ----------------------------------------------------------------
    ... 1. GET CONTOURS
        %-----------------------------------------------------------------
    
    % TURN FRAME INTO SURFACE
    for x = 1:size(frame,1)
        for y = 1:size(frame,2)
            xs(x,y) = x;
            ys(x,y) = y;
            zs(x,y) = frame(x,y);
        end
    end
    
    % INTERSECT WITH PLANE
    [M,~] = contourf( xs,ys,zs, [thresh thresh]);
    close;
    
    %TURN INTO COORDS
    if ~isempty(M)
        [x,y,z] = C2xyz(M) ; %get coordinates of the intersection
        ncoor{fri}= length(x);
        Xcoor{fri} = x;
        Ycoor{fri} = y;
        Zcoor{fri} = z;
        
    else
        ncoor{fri}= 0;
        Xcoor{fri} = 0;
        Ycoor{fri} = 0;
        Zcoor{fri} = 0;
    end
    clear x y z M
    
    % ----------------------------------------------------------------
    ... 2. GET AREA INSIDE ROIS
        %-----------------------------------------------------------------
    
    active_mask = frame> thresh;
    ti = ti+1;
    pix(ti).ms = VSDI.timebase(frameidx(fri));
    
    for roi = roiA_idx
        roimask = squeeze(roiA_mask(:,:,roi));
        
        if highdefinition
            roimask = interp2(single(roimask),5, 'cubic');
        end
        
        intersection = roimask & active_mask;
        
        roiname = roiA_labels(roi);
        roiname = roiname{1};
        pix(ti).(roiname) = sum(intersection(:));
        
        %                 % DEBUG CONTROL ······································
        %                 figure;
        %                 subplot(3,1,1)
        %                 ax0 = axes;
        %                 imagesc(back); colormap(ax0, bone)
        %                 title(['t' num2str(pix(ti).ms) 'ms' num2str(sum(roimask(:))) 'roi pix ; ' num2str(sum(active_mask(:))) 'active pix; ' num2str(sum(intersection(:))) 'intersected pix'], 'Interpreter', 'none')
        %
        %                 ax1 = axes;
        %                 imagesc(roimask,'alphadata',roimask*0.5); colormap(ax1,parula)
        %                 ax1.Visible = 'off';
        %                 linkprop([ax0 ax1],'Position');
        %
        %                 ax2 = axes;
        %                 imagesc(active_mask,'alphadata',active_mask*0.5); colormap(ax2,jet)
        %                 ax2.Visible = 'off';
        %                 linkprop([ax0 ax2],'Position');
        %
        %                 ax3 = axes;
        %                 imagesc(intersection,'alphadata',intersection*0.5); colormap(ax3,flipud(jet))
        %                 ax3.Visible = 'off';
        %                 linkprop([ax0 ax3],'Position');
        %                 % END OF DEBUG CONTROL ······································
        
    end % for roi
    
end % for fri




% ----------------------------------------------------------------
... PLOT AREA
    % ------------------------------------------------------------

% DIVIDE IN R vs L from each roi
idx2 = []; idx4 = [];
for roi = roiA_idx
    roiname = roiA_labels(roi);
    roiname = roiname{1};
    % rois with 'R'
    if strfind(roiname, '4')
        idx4 = [idx4  roi];
    elseif strfind(roiname, '2')
        idx2 = [idx2  roi];
    end
end % for roi

% Get labels for x axis (add an extra number for each extreme)
ms_labels(1) = NaN;
ms_labels(2:nframes+1) = [pix.ms];
ms_labels(end+1) = NaN;

% Dm4
if save_timemap
    lineW = 3;
else
    lineW = 1.3;
end
figure
subplot(2,2,1)
plot(1:nframes, [pix.dm4_R], '-o', 'linewidth', lineW); hold on;
plot(1:nframes, [pix.dm4_L], '-o', 'linewidth', lineW);
legend ({'contra'; 'ipsi'}, 'Location', 'Southeast')
xticks(1:nframes)
xticklabels(ms_labels(2:end-1))
xtickangle(90);
xlabel('ms')
ylabel ('nº pixels ')
title('pixels above threshold in Dm4')

% Draw Dm4 rois considered
ax01= subplot(2,2,2);
imagesc(back); colormap(ax01, bone)
axis image
% title('roi area')
for roi = idx4
    roimask = roiA_mask(:,:,roi);
    ax(roi) = axes;
    imagesc(roimask,'alphadata',roimask*0.5); colormap(ax(roi),parula)
    axis image
    ax(roi).Visible = 'off';
    linkprop([ax01 ax(roi)],'Position');
end % for roi
ax01.Visible = 'off';

% Dm2
subplot(2,2,3)
plot(1:nframes, [pix.dm2_R], '-o', 'linewidth', lineW); hold on;
plot(1:nframes, [pix.dm2_L], '-o', 'linewidth', lineW);
legend ({'contra'; 'ipsi'}, 'Location', 'Southeast')
xlim([1 nframes+1])
xticks(1:nframes)
xticklabels(ms_labels(2:end-1))
xtickangle(90);
xlabel('ms')
ylabel ('nº pixels ')
title('pixels above threshold in Dm2')

% Draw Dm2 rois considered
ax02 = subplot(2,2,4);
imagesc(back); colormap(ax02, bone)
axis image

% title('roi area')
for roi = idx2
    roimask = roiA_mask(:,:,roi);
    ax(roi) = axes;
    imagesc(roimask,'alphadata',roimask*0.5); colormap(ax(roi),parula)
    axis image
    ax(roi).Visible = 'off';
    linkprop([ax02 ax(roi)],'Position');
end % for roi
ax02.Visible = 'off';

savename= ['TIMEMAP' num2str(VSDI.ref) ref_movie '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'thresh' num2str(fact_thresh) 'perc-'  'pre' num2str(fr_pre) 'post' num2str(fr_post) 'fr' num2str(nframes) timespan_mode '-' num2str(fact_thresh*100) '%' thresh_mode];


if save_timemap
    save_currentfig_hi(savein, [savename, '1'] , 300, 'jpg')
    %         set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
end

% %% ----------------------------------------------------------------
% ... PLOT ROI TO CHECK
%     % ------------------------------------------------------------
% figure;
% ax0 = axes;
% imagesc(back); colormap(ax0, bone)
% title('roi area')
% for roi = roiA_idx
%     roimask = roiA_mask(:,:,roi);
%     ax(roi) = axes;
%     imagesc(roimask,'alphadata',roimask*0.5); colormap(ax(roi),parula)
%     ax(roi).Visible = 'off';
%     linkprop([ax0 ax(roi)],'Position');
%
% end % for roi

% ----------------------------------------------------------------
... PLOT CONTOURS
    %------------------------------------------------------------
figure
imagesc(backgr); colormap bone
axis image
hold on
timecmap = flipud(jet(numel(frameidx)));
for fri = 1:nframes
    
    if ncoor{fri}> 0
        for n= 1:ncoor{fri}
            x = Xcoor{fri}{n};
            y= Ycoor{fri}{n};
            
            plot(y, x, 'linewidth' , 2.5, 'color', timecmap(fri,:))
        end
    end
    clear x y z
    
end % for fri
mA = num2str(round(VSDI.condition(sel_trials(1),4),2));
titulo = [num2str(VSDI.ref) '.Cond' num2str(condi) '(' mA 'mA).' 'thresh' num2str(fact_thresh) 'perc-'  'pre' num2str(fr_pre) 'post' num2str(fr_post) 'fr' num2str(nframes) timespan_mode '-' num2str(fact_thresh*100) '%' thresh_mode];
sgtitle(titulo, 'Interpreter' , 'none')

%----------------------------------------------------------------
... FAKE PLOT TO GET TIME LEGEND
    %------------------------------------------------------------

timelabels = num2str(VSDI.timebase(frameidx));

for i = 1:nframes
    if ncoor{fri}> 0
        hi = plot([1 1], 'Color', timecmap(i,:), 'linewidth', 1.3);
        h(i) = hi(1);
    end
            if i == 1; hold on; end

end
lgnd= legend;
set(lgnd,'edgecolor', 'none', 'location', 'northeast', 'string', timelabels); %set(lgnd,'color','none', 'textcolor', 'w')
%----------------------------------------------------------------
... SAVE FIGURES
    %------------------------------------------------------------
% ATT:  automatical saving does not respect the legend configuration
% (save definite images manually)


if save_timemap
    save_currentfig_hi(savein, [savename, '2'], 300, 'jpg')
    %         set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
    close
end

blob()


%% Update history:
% 12/030/22 Created
% 15/03/22 - fixed bug (center in peak, method3)