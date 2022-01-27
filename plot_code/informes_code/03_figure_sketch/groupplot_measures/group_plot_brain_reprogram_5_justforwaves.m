% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% BLANK SUBSTRACTION OF THE WAVES

% average movie
% extract wave for each pixel and make the measure maps (for all
% conditions): for condition for xi for yi

% z-spatial all conditions for each fish
% extract roi from z-spatially maps >>> definite measures

clear
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new.mat') % ATTENTION: select cases manually (there are repetitions)

% ///////////////////////////////////////////////////////////
% SETTINGS

sourcedata = 'normal';
% sourcedata = '%F';
% sourcedata = 'blank-s'; % BLANK-SUBSTRACTION
% sourcedata = '%deltaF blank-s';

% selroinames = {'dm4m_R',  'dm2_R'};
% selroinames = {'dm4m_R2',  'dm2_R2'};

% selroinames = {'dm4m_R',  'dm2_R', 'dm3_R' ,'dm1_R','dldm_R'};%dm3,
selroinames = {'dm4m_R2',  'dm2_R2', 'dm3_R2' ,'dm1_R','dldm_R2'};%dm3,

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_17filt5' ;
% ref_movie= '_12filt5' ;

activ_unit = 'diffF'; % @ MANUALLY SET
analysisref = 'group2_'; %extra info for the name

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/Zscore' ;%@ SET

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; % DO NOT TOUCH: ms Range of analysis

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 600]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-25 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350];

% feedf.noise.fr_abovenoise = 30;
% feedf.noise.SDfactor = 2;
% feedf.noise.basel = [-200 0]; %it won't be used by the function, but will be used to manually get the baseline
feedf.method = 'movsum';

% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{4,1} = 'wmean (ms)';
params{4,2} = feedf.window.wmean;
% END

% ///////////////////////////////////////////////////////////////

% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');

addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));
% END USER_SETTINGS

 reject_on = [1] ;  %@ SET
    
    %----------------------------------------------------------------
    % @SET: REJECT SETTINGS
    %----------------------------------------------------------------
    
    % Subsettings:
    setting.manual_reject = 1; %@ SET
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET+
    setting.force_include = 0; %@ SET
    
    
    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    i = 2; % counter for long format rows: the first will be the labels
    si = 1; % counter for subjects list
    
    for suji  = 2%  [1:4 9] %1:size(groupplot,1)
        
        nfish = groupplot{suji,1};
        
        VSDI = TORus('load', nfish);
        VSDmov = TORus('loadmovie',nfish,ref_movie);
        movies = VSDmov.data ;
        F0 = VSDmov.F0;
%         cond_lohi = groupplot{suji,3}; 
        cond_lohi = [401 402 403]; %CHOOSE ONLY 3

%         temp = VSDI.condition(:,1);
%         temp(any(isnan(temp), 2), :) = []; %to delete rows with NaN
%         cond_lohi = unique(temp, 'rows');

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
        
        % -------------------------------------------
        % FIRST LOOP THROUGH CONDITIONS AND PIXELS TO GET PIXEL-MEASURES (of all regions, i.e., a
        % common threshold to all)
        % -------------------------------------------
        % we get and store the value in this first loop to get the max
        % value, and in the second loop we use that max-val as
        % threshold
        
        ci = 0;
        for condition =  makeRow(cond_lohi)
            ci = ci+1;
            
            cond_blank = force0ending(condition);
            % -------------------------------------------
            % SELECT CASES  AND AVERAGE MOVIE
            % -------------------------------------------
            sel_trials = find(VSDI.condition(:,1)==condition);
            sel_blank = find(VSDI.condition(:,1)==cond_blank);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
                sel_blank = setdiff(sel_blank, rejectidx);
                
            end
            
            % --------------------------------------------------------------------------
            % GET AVERAGE MOVIE. Use timerange set
            % --------------------------------------------------------------------------
            
            switch sourcedata
                
                case 'normal'
                    refcase = '';
                    
                    movieave = mean(movies(:,:,idxrange,sel_trials),4);
                    movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    ...(that normally corresponds to the background)
                        
                
                case  'blank-s'
                    refcase = 'blankS';
                    
                    movieave = mean(movies(:,:,idxrange,sel_trials),4);
                    movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    
                    movieblank = mean(movies(:,:,idxrange,sel_blank),4);
                    movieblank(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    %                         blankwave = roi_TSave(movieblank,roimask);
                    %                         roiwave =  roi_TSave(movieave,roimask)-blankwave;
                    %             case '%deltaF blank-s'
                    %                 refcase = '%deltaF blanks';
                    %                         maxval = max(roiwave(:));
                    
            end %case
            
         %----------------------------------------------------------------
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
        
            
            % -------------------------------------------------------
            % CALCULATE WAVE FOR EACH ROI 
            % -------------------------------------------------------
            for roii = makeRow(selroi)
                roimask = masks(:,:,roii);
%                 allroi_waves(:,roii,ci) = roi_TSave(movieave,roimask);
                meanF0 = squeeze(mean(F0(:,:,sel_trials),3));
                allroi_waves(:,roii,ci) = roi_TSave_percF_roiwise(movieave,roimask, meanF0);
            end %for roi_i

            
            
        end %for condition
        
        
        %% PLOT MAPS
        savename = [num2str(VSDI.ref) '_' ref_movie '_' refcase '_rej' num2str(reject_on) '_'];
%         
        ncond = length(cond_lohi);

        % ------------------------------------------------------------------
        % PLOT ALL ROI
        % ------------------------------------------------------------------
        ci = 0;
        for condition =  makeRow(cond_lohi)
            ci = ci+1;
            subplot(2,ncond+1,ci+1);
                plot(timebase_adj , allroi_waves(:,selroi,ci), 'linewidth', 2);
                hold on
                display(condition)
                ylabel('%F')
                legend(selroinames{:}, 'Location', 'southeast')
                title(num2str(condition))
%                 idxcond = find(VSDI.condition(:,1) == condition,1, 'first');
%                 condlabel{ci}= [num2str(VSDI.condition(idxcond,4)) 'mA'];        ax1 = subplot(1,4,3);
ax1 =subplot(1,ncond+1,1);

        switch roikind
            case 'circle'
                centers = VSDI.roi.circle.center(selroi, :) ;
                roicirc_preview_multiple(VSDI.crop.preview, centers, VSDI.roi.circle.R, ax1);
                
            case 'anat'
                roi_preview_multiple(VSDI.crop.preview, VSDI.roi.manual_poly(selroi,:), ax1);
        end
                sgtitle([num2str(VSDI.ref), ref_movie, 'rej', num2str(reject_on)])
        end % for condition
                end      

        localname = [num2str(VSDI.ref), ref_movie,'_ROIWAVES_',roikind,selroinames{:}, ref_movie ];
%         saveas(gcf, fullfile(savein, [localname '.jpg']), 'jpg')
%         close
        
        end

%         end