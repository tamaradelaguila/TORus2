clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
nfish = 8;
cond_codes = [2001:2003]; % DO NOT INCLUDE THE BLANK CONDITION (or the threshold will fail)

timespan_mode =  'fixed_timerange' ;  % 'peak' 'lowerthresh' 'fixed_timerange'
manual_timespan = [130 300]; % if fixed_timerange

save_timemap = 1;
highdefinition = 0;

% FOT TILES AND EARLY-PEAK FRAMES
fact_thresh =0.25; % @SET : limits parameters

movie_ref = '_18filt6'; % input movie '_17filt5'
% movie_ref = '_17filt5'; % input movie '_17filt5'

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/def_figs/timemaps/'; % CAHNGEEEEEE

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on= 3;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.force_include = 0; %

%----------------------------------------------------------------
% @SET: for roi-waves plot
%----------------------------------------------------------------

refroiname = 'dm4m_R2'; % respect which the start/end time of the tiles will be set

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
rej = 'reject' ;
if reject_on > 1
    rej = [rej num2str(reject_on)];
end

rejectidx = [];

if setting.manual_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
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
% SELECT ROI
%----------------------------------------------------------------
switch roikind
    case 'circle'
        selroi =name2idx(refroiname, VSDI.roi.labels_circ);
        roilabels = VSDI.roi.labels_circ;
        masks =  VSDI.roi.circle.mask;
        
    case 'anat'
        selroi =name2idx(refroiname, VSDI.roi.labels);
        roilabels = VSDI.roi.labels;
        masks = VSDI.roi.manual_mask;
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


%-------------------------------------------------------------------
% GET MAX FROM ALL CONDITIONS TO SET AS REPRESENTATION'S THRESHOLD
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
    movieave = mean(VSDmov.data(:,:,:,sel_trials),4);
    movieave(:,:,end) = back; %clean non-blured background
    
    % SMOOTH DATA TO GET THE MAX VALUE
    tempmax = movmean(movieave(:,:,1:end-1),5); %temporal smooth to get the max value par: 3
    
    maxall(ci) = max(tempmax(:));
end

maxval = max(maxall); % ATT !!!

%----------------------------------------------------------------
... LOOP TO GET TILES, EARLY-PEAK FRAMES AND WAVES
%----------------------------------------------------------------

for condi = makeRow(cond_codes)
    
    %----------------------------------------------------------------
    ... GET AVERAGE MOVIE
        %----------------------------------------------------------------
    [sel_trials] = find(VSDI.condition(:,1)==condi);
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    
    back = VSDI.backgr(:,:,VSDI.nonanidx(1));
    
    %to plot single trial
    movieave = mean(VSDmov.data(:,:,:,sel_trials),4);
    movieave(:,:,end) = back; %clean non-blured background
    
    
    %----------------------------------------------------------------
    % TILES
    %----------------------------------------------------------------
    
    
%     tileset.clims = [-maxval*fact_clim maxval*fact_clim];
    thresh = [-maxval*fact_thresh maxval*fact_thresh];
    
    % -------------------------------------------------------
    % CALCULATE dF WAVE FOR dm4 TO ESTABLISH START-END TIMES
    % -------------------------------------------------------
    
    idxDm4 =name2idx(refroiname, roilabels);
    roimask = masks(:,:,idxDm4);
    
    dm4_wave = roi_TSave(movieave,roimask);
    dm4_wave = movmean(dm4_wave ,5);
    wave_timebase = VSDI.timebase;
    
    fr_pre = 0; %nº of frames before rising the threshold (to set the initial time)
    fr_post = 2;

    idx.start = find(dm4_wave > thresh(2), 1, 'first')- fr_pre;
    if isempty (idx.start)
        error('The threshold is too high (lower the threshold by changing the value of "fact_thresh"')
    end

    
    switch timespan_mode
        case 'lowerthresh'
            idx.end= find(dm4_wave > thresh(2), 1, 'last') + fr_post;
        case 'peak'
            temp = devo_peak2peak(dm4_wave, wave_timebase, feedf.window, [],  'movsum' , 0, 0);
            idx.end =temp.peakidx(2);
        case   'fixed_timerange'
            idx.start = dsearchn(VSDI.timebase, manual_timespan(1));
            idx.end = dsearchn(VSDI.timebase, manual_timespan(end));
    end
    
    if isempty (idx.start)
        error('The threshold is too high (lower the threshold by changing the value of "fact_thresh")')
    end
    
%     titulo = [num2str(VSDI.ref) 'cond' num2str(condi) '.clim=' num2str(round(tileset.clims(2),1)) '(' num2str(fact_clim) '%)', '.thresh=' num2str(round(thresh(2),1)) '(' num2str(fact_thresh) '%)'];
%     sgtitle(titulo)
    
    
    % /////////////////////////////////////////////////////////////////
    %----------------------------------------------------------------
    ... PLOT EARLY-PEAK
        %------------------------------------------------------------
    % /////////////////////////////////////////////////////////////////
    nframes = 6;
    frameidx = round( linspace(idx.start,idx.end, nframes));

    backgr = movieave(:,:,end);
        if highdefinition
        backgr = interp2(backgr,5, 'cubic');
        end

    % FUTURE FUNCTION
    for fri = 1:numel(frameidx)

        frame = movieave(:,:,frameidx(fri));
        
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
        
        %             figure
        %             s= surf(xs,ys,frame); colormap bone
        %             h= surf(xs,ys,frame, 'FaceAlpha',0.5)
        %             s.EdgeColor = 'none';
        %             view([-37.5 -30])
        %             figure
        % %             s= surf(xs,ys,zs); colormap bone
        %             h= surf(xs,ys,zs, 'FaceAlpha',0.5);
        %             s.EdgeColor = 'none';
        %             view([-37.5 -30])
        %
        
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
        %             % PLOT ON THE IMAGE
        %             figure
        %             imagesc(backgr); colormap bone
        %             hold on
        %
        %             for ii = 1:L
        %                 xx= x{ii};
        %                 yy = y{ii};
        %                 plot(yy, xx, 'linewidth' , 1.5)
        %             end
        %
    end % for fri
    
    %----------------------------------------------------------------
    ... PLOT CONTOURS
        %------------------------------------------------------------
    figure
    imagesc(backgr); colormap bone
    axis image
    hold on
    timecmap = flipud(jet(numel(frameidx)));
    for fri = 1:numel(frameidx)
        
        if ncoor{fri}> 0
            for n= 1:ncoor{fri}
                x = Xcoor{fri}{n};
                y= Ycoor{fri}{n};
                
                plot(y, x, 'linewidth' , 1.5, 'color', timecmap(fri,:))
            end
        end
        clear x y z
        
    end % for fri
    mA = num2str(round(VSDI.condition(sel_trials(1),4),2));
    titulo = [num2str(VSDI.ref) '.Cond' num2str(condi) '(' mA 'mA).'  'thresh=' num2str(round(thresh(2),1)) '(' num2str(fact_thresh*100) '%)'];
    sgtitle(titulo)
    
    %----------------------------------------------------------------
    ... FAKE PLOT TO GET TIME LEGEND
        %------------------------------------------------------------
    
    timelabels = [num2str(VSDI.timebase(frameidx))];
    clear h
    for i = 1:numel(frameidx)
        hi = plot([1 1], 'Color', timecmap(i,:), 'linewidth', 1.3);
        h(i) = hi(1);
%         if i == 1; hold on; end
    end
    
    lgnd = legend(h, timelabels);
    set(lgnd,'edgecolor', 'none', 'location', 'northeast',  'color', 'none', 'textcolor', 'w'); %set(lgnd,'color','none')
    
    %----------------------------------------------------------------
    ... FAKE PLOT TO GET TIME LEGEND
    %------------------------------------------------------------
    % ATT:  automatical saving does not respect the legend configuration
    % (save definite images manually)
    
    savename= ['TIMEMAP' num2str(VSDI.ref) movie_ref '_cond' num2str(condi) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'thresh' num2str(fact_thresh) 'perc-'  'pre' num2str(fr_pre) 'post' num2str(fr_post) 'fr' num2str(nframes) '.jpg'];
    
    if save_timemap
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
        print(fullfile(savein, [savename]),'-r600','-djpeg') % prints it as you see them
    end
    
end % for condi
blob()
