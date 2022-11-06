% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% BLANK SUBSTRACTION OF THE WAVES

% average movie
% extract wave for each pixel and make the measure maps (for all
% conditions): for condition for xi for yi

% z-spatial all conditions for each fish
% extract roi from z-spatially maps >>> definite measures
close all
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus';
user_settings
cd(W)

% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new.mat') % ATTENTION: select cases manually (there are repetitions)

% FOR NEW BOXPLOTS (no-pain vs pain)
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot4_boxplot.mat') % ATTENTION: select cases manually (there are repetitions)

% ...FOR 'EFFECT OF INTENSITY' BOXPLOT:
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot5_boxplot_Ieffect.mat') % ATTENTION: select cases manually (there are repetitions)

% ...FOR 'CHORRITO' BOXPLOT:
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot6_chorrito.mat') % ATTENTION: select cases manually (there are repetitions)

% ///////////////////////////////////////////////////////////
% SETTINGS
selroinames = {'dm4m_R2',  'dm2_R2'};

% Reference for the code
latency_refroi = 'dm2_R2';% to get percentage for threshold computation
latency_thresh = 0.8; % 80% of the refroi value

roikind = 'circle'; %
% roikind = 'anat';

dataunits = 'dF'; % '%F' 'dF'
switch dataunits
    case 'dF'
        movie_ref = '_18filt6'; % input movie '_17filt5'
    case '%F'
        movie_ref = '_21filt6';
end

ref_movie= '_18filt6';% '_17filt5' ; '_18filt6'
% ref_movie= '_12filt5' ;

activ_unit = 'dF'; % @ MANUALLY SET (just for display purposes)
analysisref = 'new4_group19_dolor_n5'; % new4_group16_chorrito_n3  new4_group19_dolor_n5
reject_on = 3; 

getR = 1;

savewaves = 1;
% savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/Zscore' ;%@ SET
% savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/04_figure_sketch/groupplot_measures/Zscore' ;%@ SET
% savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/04_figure_sketch2/groupplot_measures/'; %10/09/22
savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit';

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 1200]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-25 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350];

% feedf.noise.fr_abovenoise = 30;
% feedf.noise.SDfactor = 2;
% feedf.noise.basel = [-200 0]; %it won't be used by the function, but will be used to manually get the baseline
feedf.method = 'movsum';

% flagslope = 1;

% Params
slope.window = [0, 200]; %ms

% END USER_SETTINGS
% ///////////////////////////////////////////////////////////////

%% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{5,1} = 'slope win (ms)';
params{5,2} = slope.window;

params{6,1} = ['onset lat (ms) ' num2str(latency_thresh*100) '% from max of:' latency_refroi];

params{8,1} = date;
params{8,3} = ['source:'  mfilename('fullpath')];


%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

% Subsettings:
setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 0; %@ SET+
setting.forcein = 0; %@ SET

% SUBJECTS SELECTION (ROWSA FROM THE 'groupplot' CELL STRUCTURE
% FOR NO-PAIN
...................
    
% FOR PAIN
...................
    if strcmpi( analysisref, 'new4_group17_pulsito') % del panel original: para dm2+dm4 vs blank
    sel_subjects = [24 25 26 29 30];
    elseif strcmpi( analysisref, 'new4_group19_dolor_n5') % del panel original: para dm2+dm4 vs blank
        sel_subjects = [31:35];
    elseif strcmpi( analysisref, 'new4_group19_dolor_n4') % del panel original: para dm2+dm4 vs blank
        sel_subjects = [31:34];
        
        % FOR INTENSITY EFFECT BOXPLOTS
        ...................
            % from: load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot5_boxplot_Ieffect.mat') % ATTENTION: select cases manually (there are repetitions)
        
    elseif strcmpi( analysisref, 'new4_group20_3intens') % del panel original: para dm2+dm4 vs blank (n=6)
        sel_subjects = [11:16]; %n=4
        
        % FOR CHORRITO
        % ...................
        % from: load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/04_figure_sketch/groupplot_measures/groupplot6_chorrito.mat') % ATTENTION: select cases manually (there are repetitions)
    elseif strcmpi( analysisref, 'new4_group16_chorrito_n3') % del panel original: para dm2+dm4 vs blank
        sel_subjects = 1:3;
    elseif strcmpi( analysisref, 'new4_group16_chorrito_n4') % del panel original: para dm2+dm4 vs blank
        sel_subjects = 1:4;
    end
    
    
    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    i = 2; % counter for long format rows: the first will be the labels
    si = 0; % counter for subjects list
    
    for suji  =  sel_subjects %1:size(groupplot,1) SELECT included fish+condition
        
        si = si+1;
        nfish = groupplot{suji,1};
        cond_list = groupplot{suji,3};
        cond_list = cond_list(2:end);
        
        VSDI = TORus('load', nfish);
        VSDmov = TORus('loadmovie',nfish,ref_movie);
        movies = VSDmov.data ;
        F0 = VSDmov.F0;
        
        %----------------------------------------------------------------
        % CONTROL ROI PICTURE
        %----------------------------------------------------------------
        %     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R);
        %     title([num2str(VSDI.ref) 'roi preview:' selroinames{:}])
        %     saveas(gcf, [num2str(VSDI.ref)'roipreview'], 'jpg')
        
        
        %----------------------------------------------------------------
        % GET INDEXES OF TIMERANGE
        %----------------------------------------------------------------
        idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
        idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values
        
        timebase_adj = VSDI.timebase(idxrange);
        
        % get indexes for slope
        %         slope.windowidx = dsearchn(timebase_adj, makeCol(slope.window));
        %         slope.windowidx = [slope.windowidx(1) slope.windowidx(end)];
        
        %----------------------------------------------------------------
        % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
        %----------------------------------------------------------------
        
        rejectidx = [];
        rejectidx  = compute_rejectidx(VSDI, reject_on, setting);
        
        %         params{8,1} = 'rejected';
        
        %         params{8,1} = rejectidx';
        
        % Sampling frequency for lowpass flter
        %----------------------------------------------------------------
        
        fs = 1000/(VSDI.info.stime);
        
        
        %% ----------------------------------------------------------------
        % SELECT ROI
        %----------------------------------------------------------------
        
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
        
        % REFERENCE ROI
        refidx=name2idx(latency_refroi, roilabels);
        ref_roimask = masks(:,:,refidx);
        
        
        
        % -------------------------------------------
        % FIRST LOOP THROUGH CONDITIONS AND PIXELS TO GET PIXEL-MEASURES
        % -------------------------------------------
        
        ci = 0;
        for condition =  makeRow(cond_list)
            ci = ci+1;
            
            %             cond_blank = 0; % from june2022 on... when there is a single blank trial
            
            % -------------------------------------------
            % SELECT CASES  AND AVERAGE MOVIE
            % -------------------------------------------
            sel_trials = find(VSDI.condition(:,1)==condition);
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            % --------------------------------------------------------------------------
            % GET AVERAGE MOVIE. Use timerange set
            % --------------------------------------------------------------------------
            
            movieave = mean(movies(:,:,idxrange,sel_trials),4);
            movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
            ...(that normally corresponds to the background)
                
        % --------------------------------------------------------------------------
        % GET AVERAGE MOVIE. Use timerange set
        % --------------------------------------------------------------------------
        % Refroi wave
        reflat_wave = roi_TSave(movieave,ref_roimask);
        %                 reflat_wave = movmean(reflat_wave ,5);
        reflat_wave = lowpass(reflat_wave,10,fs);
        
        %Get peak of that wave, multiply by factor and set as threshold
        %for latency calculation
        
        % method 1 - peak finding
                    temp0 = devo_peak2peak(reflat_wave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
        %             cond_thresh = temp0.peakminusbasel*latency_thresh;
        
        % method 2:
        [~,idxpeak2] = max(movsum(reflat_wave, 15));
        cond_thresh2 = reflat_wave(idxpeak2) *latency_thresh;
        feedf.window.risingthres  =  cond_thresh2;
        
        % method 3:
        idxpeak = temp0.peakidx(2);
        cond_thresh = reflat_wave(idxpeak) *latency_thresh;
        
        % TEST
        % .................................................................
        figure
        plot(reflat_wave, 'LineWidth', 1.8)
        xline(idxpeak2, 'r-', {'max movsum'}); xline(idxpeak, 'g-', {'devo_peak2peak'});
        yline(cond_thresh2, 'r--'); yline(cond_thresh, 'g--');
        title('two methods for threshold determination')
        saveas(gcf, fullfile(savein, ['TEST_PEAK_' num2str(VSDI.ref) '_cond' num2str(condition)  '.jpg']), 'jpg')
        close
        
        %% ---------------------------------------------------------------
        ... LATENCY FROM WAVES:
            % ---------------------------------------------------------------
        % ploti = 0;
        
        % Write to
        coli = 1; %counter for columns
        fishwaves(:,coli) = timebase_adj;
        labels{1,coli} =   'ms';
        
        figure
        for roi_i = makeRow(selroi)
            %     ploti = ploti+1;
            %     subplot(1,length(selroi), ploti)
            roimask= masks(:,:,roi_i);
            
            roiwave =  roi_TSave(movieave,roimask);
            roiwave = lowpass(roiwave, 10, fs);
            
            %     temp = devo_peak2peak(roiwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
            idxlat = find(roiwave> cond_thresh, 1, 'first');
            peaklat(roi_i,ci) = timebase_adj(idxlat);
            
            %     % PLOT TO CHECK
            plot(timebase_adj, roiwave, 'linewidth', 1.8, 'displayname',  roilabels{roi_i})
            xline(peaklat(roi_i,ci), '-', {['lat=' num2str(round(peaklat(roi_i,ci))), 's']}, 'handlevisibility', 'off');
            yline(cond_thresh, '-', {'threshold'}, 'handlevisibility', 'off') ;
            
            hold on
            ylabel('dF'); xlabel('ms')
            
            
            clear temp
            
            % Get matrix to later write in excel
            coli = coli+1;
            fishwaves(:,coli) = roiwave;
            fishwaves(1,coli+numel(selroi)) = peaklat(roi_i,ci); % print latency
            labels{1,coli} =   roilabels{roi_i};
            labels{1,coli+numel(selroi)} =  ['lat' roilabels{roi_i}];
            
        end
        legend
        hold off
        
        sgtitle([num2str(VSDI.ref), movie_ref, 'rej', num2str(reject_on), '-cond:' num2str(condition)])
        
        savename= ['LATENCY_WAVES' num2str(VSDI.ref) movie_ref '_cond' num2str(condition) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' ];
        
        if savewaves
            saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
            print(fullfile(savein, [savename '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them %STILL TO TEST!
            
            close
        end
        
        
        
        % sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Peaklatency from wave .Cond:' num2str(condition)])
        %
        % localname = ['Peaklat_from_avewave_cond' num2str(condition)] ;
        % savename = [num2str(VSDI.ref) '_' ref_movie '_rej' num2str(reject_on) '_' ];
        %
        % saveas(gcf, fullfile(savein, [savename localname '.jpg']), 'jpg') % prints it as you see them
        % close
        
        clear cond_thresh
        end %for condition
        
        
        %% -------------------------------------------
        % EXPORT FOR R
        % --------------------------------------------
        % one excel for all group, one sheet for each fish/condition
        if getR
            excelname = fullfile(savein, [analysisref '_RAWWAVES_latency' ref_movie '_rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls']);
            
            for col = 1:numel(labels)
                labels{1,col}= labels{col};
            end
            
            % write output (new sheet for each fish
            writematrix (fishwaves, excelname, 'sheet',  [num2str(VSDI.ref) '-' num2str(condition)])
            
            writecell (labels, excelname, 'sheet', 'labels')
            writecell (params, excelname, 'sheet','param')
            
            clear longF
        end % if getR
        
    end % for suji
    
    % LOOP TO GET NORMALIZED WAVE
    
    
    blob()
    
    % Update history:
    % 27/10/22 - FIX ISSUE WITH PEAK FINDING FOR THRESHOLD: should not be
    % the peakminusbasel but the peak
    % 16/02/22 - Add onset-latency measure and delete 'switch refcase' (so only
    % the normal dF condition is considered). To use blank-substraction
    % methods, go to the older versions of the code
