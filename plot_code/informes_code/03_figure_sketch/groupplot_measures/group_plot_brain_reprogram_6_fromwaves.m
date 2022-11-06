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

% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_latency.mat')
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new3.mat') % ATTENTION: select cases manually (there are repetitions)

% ///////////////////////////////////////////////////////////
% SETTINGS

% selroinames = {'dm4m_R',  'dm2_R'};
% selroinames = {'dm4m_R2',  'dm2_R2'};

% selroinames = {'dm4m_R',  'dm2_R', 'dm3_R' ,'dm1_R','dldm_R'};%dm3,
% selroinames = {'dm4m_R2',  'dm2_R2', 'dm3_R2' ,'dm1_R','dldm_R2'};%dm3,
% selroinames = {'dm4m_R',  'dm4_L', 'dm2_R' ,'dm2_L'};
% selroinames = {'dm4m_R2','dm2_R2','dm1_R','dldm_R2'}; %4roi
% selroinames = {'dm4m_R2','dm2_R2','dm1_R','dldm_R2','dm4m_L2','dm2_L2','dm1_L','dldm_L2'}; %8roi
selroinames = {'dm4m_R2', 'dm4m_L2', 'dm2_R2', 'dm2_L2'}; %4roi

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_18filt6';% '_17filt5' ; '_18filt6'
% ref_movie= '_12filt5' ;

activ_unit = 'diffF'; % @ MANUALLY SET (just for display purposes)
analysisref = 'testlateral_group8_'; % MANUALLY SET!!! extra info for the name. Set group according to the rows of 'groupplot' selected in;   for suji  =  [1 3 4 9] lateral_group7_RECOV
reject_on = 3;

onset_factor = 0.2; % of max value to set rising threshold (to onset latency)

plotwaves = 1;
savewaves = 0; % if plotwaves =1

plottest = 1; 
savetest =0;

getR = 1;
get_excel_indiv = 0; %not working at the moment

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/measures_fromwaves' ;%@ SET

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

% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{4,1} = 'wmean (ms)';
params{4,2} = feedf.window.wmean;

params{5,1} = 'slope win (ms)';
params{5,2} = slope.window;

params{6,1} = ['onset lat ' num2str(onset_factor*100) '%(ms)'];

if strcmpi( analysisref, 'panel2_group8')
%     sel_subjects = [1 2 5 6]; % from groupplot_new2.mat
        sel_subjects = [1 2 3 4]; % from groupplot_new3.mat
elseif strcmpi( analysisref, 'testlateral_group8_')
%     sel_subjects = [1 2 5 6]; % from groupplot_new2.mat
        sel_subjects = [1 2 3 4]; % from groupplot_new3.mat
        
end
    %----------------------------------------------------------------
    % @SET: REJECT SETTINGS
    %----------------------------------------------------------------
    
    % Subsettings:
    setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET+
    setting.force_include = 0; %@ SET
    

    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    i = 2; % counter for long format rows: the first will be the labels
    si = 0; % counter for subjects list
    
    for suji  =  sel_subjects %1:size(groupplot,1) SELECT included fish+condition
%         
%         si = si+1;
        nfish = groupplot{suji,1};
        condition = groupplot{suji,3};
%         
        % DEBUG !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
%         nfish = 11;  % DEBUG
%         condition = 402; 
        
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
        slope.windowidx = dsearchn(timebase_adj, makeCol(slope.window));
        slope.windowidx = [slope.windowidx(1) slope.windowidx(end)];

        %----------------------------------------------------------------
        % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
        %----------------------------------------------------------------

rejectidx = [];
reject_idx  = compute_rejectidx(VSDI, reject_on, setting);
        
        
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

        % -------------------------------------------
        % GET WAVES AND THRESHOLD
        % -------------------------------------------
                    
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
                            
            % -------------------------------------------
            % CALCULATE MEASURES for roi
            % -------------------------------------------
                    
                    for roi = makeRow(selroi)
                        roimask = masks(:,:,roi);
                            
                            roiwave = roi_TSave(movieave, roimask);
                            temp = devo_peak2peak(roiwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
                            
                            waves.peak(roi) = temp.peakminusbasel;
                            waves.wmean(roi) = temp.wmean;
                            waves.peaklat(roi)  = temp.peaklatency;
                            
                            thresh = onset_factor*waves.peak(roi);
                            feedf.window.risingthresh = thresh;
                            temp2 = devo_peak2peak(roiwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
                            temp3 = find(roiwave>onset_factor*waves.peak(roi), 1, 'first');
                            waves.onsetlat(roi) = timebase_adj(temp3);                          
                            
                            % SLOPE
                                idx0 = slope.windowidx(1);
                                idxend = slope.windowidx(end);
                                waveW = roiwave(idx0:idxend);
                                slopemean = mean(diff(waveW));
                                waves.slope(roi) = slopemean;
                                
                                
                                if plottest
                                    
                                    plot(timebase_adj,roiwave);
                                    yline(waves.peak(roi), '--' , {['peak =' num2str(waves.peak(roi))]})

                                    xline(waves.peaklat(roi), '-', {'peaklatency'})
                                    yline( onset_factor*waves.peak(roi), '--' , {[num2str(onset_factor*100) '% peak =' num2str(thresh)]})
                                    xline(waves.onsetlat(roi), '-', {['onsetLat' ]})
                                    savename = ['wavecheck_' num2str(VSDI.ref) '_' ref_movie '_rej' num2str(reject_on) '_' roilabels{roi} '.jpg'];
                                    
                                    if savewaves
                                        saveas(gcf, fullfile(savein, savename), 'jpg')
                                        close
                                    end
                                end
                             clear roiwave waveW slopemean temp temp2

                                
                    end % for roi
                    
            % get maxval for ONSET LATENCY THRESHOLD
             
               
               
               % --------------------------------------------------------------------------


        %% EXTRACT ROI MEASURES AND STORE IN LONG FORMAT (for R)
            
            if getR
                % -------------------------------------------------------
                % CALCULATE MEASURE FOR EACH ROI AND STORE IN LONG FORMAT
                % -------------------------------------------------------
                for roi_i = makeRow(selroi)
                    
                    roipeak = waves.peak(roi_i);
                    roiwmean = waves.wmean(roi_i) ;
                    roislope = waves.slope(roi_i);
                                        
                    roionsetlat =  waves.onsetlat(roi_i);
                    roipeaklat =  waves.peaklat(roi_i);

                    % FACTORS
                    longF{i,1} = suji ; %subject id
                    longF{i,2} = roi_i; % roi
                    longF{i,3} = roilabels{roi_i}; % roi
                    longF{i,4} = condition; % condition (1 = low; 2 = hi; 3 = blank)
                    % MEASURES
                    longF{i,5} = round(roipeak,2); % outputP(roi_i, condi)
                    longF{i,6} = round(roiwmean,2);%
                    longF{i,7} = round(roislope,2);%
                    longF{i,8} = round(roionsetlat);%
                    longF{i,9} = round(roipeaklat);


                    % set hemisphere factor
                    if strcmpi( roilabels{roi_i}(end), 'R') || strcmpi(roilabels{roi_i}(end-1),'R')
                        longF{i,11} = 'R';%
                    elseif strcmpi( roilabels{roi_i}(end), 'L') || strcmpi(roilabels{roi_i}(end-1),'L')
                        longF{i,11} = 'L';%
                    end
                    
                    longF{i,12} = roilabels{roi_i}(1:4);%in case we want to perform anova roi*hemisf
                    
                    i = i+1;
                    
                end %for roi_i
            end % for if getR

            
        groupplot_print(suji,:) = groupplot_latency(suji,1:end);
    end % for suji
    
    %% -------------------------------------------
    % EXPORT FOR R
    % -------------------------------------------
    if getR
        excelname = fullfile(savein, [analysisref 'Zscore_long_format_forR_group' ref_movie '_rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls']);
        labels = {'id' 'roi' 'roi n' 'cond' 'peak' 'wmean' 'slope' 'onsetlat' 'peaklat' 'hemis', 'roii'};
        
        for col = 1:numel(labels)
            longF{1,col}= labels{col};
        end
        
        % write output (new sheet for each fish
        writecell (longF, excelname, 'sheet',  [roikind ref_movie 'rej' num2str(reject_on)])
        
        writecell (labels, excelname, 'sheet', 'labels')
        writecell (groupplot_print, excelname, 'sheet','group')
        writecell (params, excelname, 'sheet','param')
        
        clear longF
    end % if getR
    
    blob()

% Update history:
% 12/03/22 - ADAPT groupplot_brain_reprogram_Zscore to extract the measures
% from the waves instead of the Pmaps
