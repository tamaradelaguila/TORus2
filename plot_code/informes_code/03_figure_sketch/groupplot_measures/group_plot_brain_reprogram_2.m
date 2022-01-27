% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% BLANK SUBSTRACTION OF THE WAVES

clear
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')

% ///////////////////////////////////////////////////////////
% SETTINGS

sourcedata = 'normal';
% sourcedata = '%deltaF';
% sourcedata = 'blank-s'; % BLANK-SUBSTRACTION
% sourcedata = '%deltaF blank-s';


% selroinames = {'dm4m_R',  'dm2_R'};
selroinames = {'dm4m_R',  'dm2_R', 'dm3_R' ,'dm1_R','dldm_R'};%dm3,

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_17filt5' ;
% ref_movie= '_12filt5' ;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures' ;%@ SET
ploton = 0;
plotbar_on = 0;
% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 600]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-25 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350];

feedf.noise.fr_abovenoise = 30;
feedf.noise.SDfactor = 2;
feedf.noise.basel = [-200 0]; %it won't be used by the function, but will be used to manually get the baseline

% Window For average-based analysis

feedf.window_ave = feedf.window;

% feedf.noise_ave = feedf.noise;
% feedf.noise_ave.SDfactor = 4;% SET differences


feedf.method = 'movsum';
% END OF FUNCTION SETTINGS

% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{4,1} = 'wmean (ms)';
params{4,2} = feedf.window.wmean;

params{5,1} = 'noise threshold';
params{5,2} = feedf.noise.SDfactor;

params{6,1} = 'SDfactor';
params{6,2} = feedf.noise.SDfactor;

params{6,1} = 'frames above';
params{6,2} = feedf.noise.fr_abovenoise;

% END


% CASE OF DATA SELECTION

switch sourcedata
    case 'normal' % normal data
        flagblank = 1;
        
    case  '%deltaF'
        flagblank = 1;

    case  'blank-s' % BLANK- SUBSTRACTED
        flagblank = 0;
        %             case '%deltaF blank-s'
        %                 refcase = '%deltaF blanks';
end

% MATCH NUMBER OF CONDITIONS WITH THE NEED TO CALCULATE THE BLANK CONDITION to consider (the blank-trial needs to be
% stored in the last position in the corresponding cell in the
% 'groupplot'cell (groupplot.mat)

if flagblank
    n = 3;
else
    n = 2; % when blank-substracted, only the two first (non-blank) conditions from the list will be computed
end


% ///////////////////////////////////////////////////////////////

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

for reject_on = [0 1]  %@ SET
    % Subsettings:
    setting.manual_reject = 0; %@ SET
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET+
    setting.force_include = 0; %@ SET
    
    
    % USER SETTINGS
    path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
    tempsep.idcs   = strfind(path.rootpath,'/');
    tempsep.newdir = path.rootpath(1:tempsep.idcs(end));
    
    path.data = fullfile(path.rootpath, 'data');
    path.grouplist = fullfile(path.rootpath);
    path.list =fullfile(path.rootpath, 'data','BVlists');
    
    addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));
    % END USER_SETTINGS
    
    % BUILD MEASURES GROUP -LONGFORMAT- MATRIX
    i = 2; j =2;
    
    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    
    for suji  =  1:size(groupplot,1)
        
        nfish = groupplot{suji,1};
        cond_lohi = groupplot{suji,3};
        
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
        % SELECT ROI
        %----------------------------------------------------------------
        
        switch roikind
            case 'circle'
                selroi =name2idx(selroinames, VSDI.roi.labels_circ);
                roilabels = VSDI.roi.labels_circ;
            case 'anat'
                selroi =name2idx(selroinames, VSDI.roi.labels);
                roilabels = VSDI.roi.labels;
        end
        
        %----------------------------------------------------------------
        % GET INDEXES OF TIMERANGE
        %----------------------------------------------------------------
        idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
        idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values
        
        % ... AND ALSO OF BASELINE FOR NOISE-THRESHOLD DETERMINATION
        baselrange =  dsearchn(makeCol(VSDI.timebase), makeCol(feedf.noise.basel));
        baselrange = baselrange(1):baselrange(end);
        
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
        % FIRST LOOP THROUGH CONDITIONS AND ROIS TO GET MEASURES AND LATENCY RESPECT TO THE MAX PEAK (of all regions, i.e., a
        % common threshold to all)
        % -------------------------------------------
        % we get and store the value in this first loop to get the max
        % value, and in the second loop we use that max-val as
        % threshold
        
        thesh = [];
        for condi =   1:n
            
            condition = cond_lohi(condi);
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
            
            movieave = mean(movies(:,:,idxrange,sel_trials),4);
            
            movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
            ...(that normally corresponds to the background)
                
        backave = mean(movies(:,:,idxrange,sel_trials),4);
        
        for roi_i =  makeRow(selroi)
            
            % STORE THE WAVES TO AVOID RECOMPUTE IT IN THE NEXT LOOP
            
            switch roikind
                case 'circle'
                    roimask = VSDI.roi.circle.mask(:,:,roi_i);
                    
                    %             anamask = VSDI.roi.manual_mask(:,:,1);
                case 'anat'
                    roimask = VSDI.roi.manual_mask(:,:,1);
            end
            
            % --------------------------------------------------------------------------
            % STORE WAVES. Use timerange set
            
            
            switch sourcedata
                
                case 'normal'
                    roiwave =  roi_TSave(movieave,roimask);
                    refcase = '';
                
                case  '%deltaF'

                    roiwave =  roi_TSave(movieave,roimask); % GET F0!!!
                    refcase = '%F';
                    
                case  'blank-s'
                    movieblank = mean(movies(:,:,idxrange,sel_blank),4);
                    movieblank(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    blankwave = roi_TSave(movieblank,roimask);
                    roiwave =  roi_TSave(movieave,roimask)-blankwave;
                    refcase = '_blankS';
                    %             case '%deltaF blank-s'
                    %                 refcase = '%deltaF blanks';
                    
            end
            
            roiwaves(:,roi_i, condi, suji) = roiwave; % GET BACK FROM PREVIOUS LOOP
            
            
            % ----------------------------------------------------------------------------------
            % GET THE PEAK VALUE FOR AN 'ABSOLUTE THRESHOLD' DETERMINATION ('20% abs thresh)
            % ----------------------------------------------------------------------------------
            temp = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window,feedf.noise, feedf.method, 0, 0);
            thresh(roi_i, condi, suji)= temp.peakminusbasel;
            
            
        end % groupplot
        end
        
        % ----------------------------------------------------------------------------------
        % GET THE ALL-CONDITIONS BASELINE MEAN AND DEVIATION FOR (CUSTOM)
        % NOISE THRESHOLD
        % ----------------------------------------------------------------------------------
        for roi_i =  makeRow(selroi)
            
            basel = squeeze(roiwaves(baselrange,roi_i, :, suji));
            d = size(basel);
            basel = reshape( basel, [1,d(1)*d(2)]); %concatenate baselines from all conditions
            thresh_noise(roi_i,suji) = mean(basel) + std(basel)*feedf.noise.SDfactor; % CUSTOM NOISE THRESHOLD BASED ON THE BASELINE OF ALL CONDITIONS
        end
    end
    %% SECOND 'suji' LOOP: MEASURES DETERMINATION (2A) AND NORMALIZATION (2B)
    
    for suji  =  1:size(groupplot,1)
        %                     testwave = roi_TSave(movieave,roimask);
        %                     plot(testwave)
        nfish = groupplot{suji,1};
        VSDI = TORus('load', nfish); %we need it to get the info of the fish for plot title
        
        %----------------------------------------------------------------
        % SELECT ROI
        %----------------------------------------------------------------
        
        switch roikind
            case 'circle'
                selroi =name2idx(selroinames, VSDI.roi.labels_circ);
                roilabels = VSDI.roi.labels_circ;
            case 'anat'
                selroi =name2idx(selroinames, VSDI.roi.labels);
                roilabels = VSDI.roi.labels;
        end
        
        % SUB-LOOP 2A:
        
        for condi = 1:n
            
            for roi_i = makeRow(selroi)
                
                roiwave = roiwaves(:,roi_i, condi, suji) ; % GET BACK FROM PREVIOUS LOOP
                
                % -------------------------------------------
                % %dF/F0
                % -------------------------------------------
                
                % plot(roiwave(1:601))
                
                % -------------------------------------------
                % CALCULATE MEASURES
                % -------------------------------------------
                
                % peakminubaseline
                feedf.window.risingthresh = thresh_noise(roi_i, suji);
                
                temp = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window, feedf.noise, feedf.method, 0, 0);
                outputP.peakminusbaseline(roi_i, condi) = temp.peakminusbasel;
                outputP.wmean(roi_i, condi) =  temp.wmean;
                %             outputP.onsetnoise_ms(roi_i, condi) = temp.onsetnoise_ms;
                % CUSTOM NOISE THRESH
                outputP.onsetnoise_ms(roi_i, condi) =  temp.risingthreshlat;
                % initialmean
                
                
                % latency respect to the 30% of the own peak
                
                %                     % -------------------------------------------
                %                     % LATENCY 30% ONSET
                %                     % -------------------------------------------
                %                     % find the threshold
                %
                %                     output_lat(roi_i, condi) = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window,[], feedf.method, 0, 0);
                
                % -------------------------------------------
                % LATENCY RESPECT  THE 20% OF ALL PEAKS
                % -------------------------------------------
                % to find the risingtime over the absolute threshold
                local = thresh(:,:,suji); %  values of all conditions and roi_i for the subject...
                feedf.window.risingthresh =  max(local(:))*0.2; %...to set the function feed parameter (absolute threshold)
                % specific inputs to use an absolute threshold
                output_abs(roi_i, condi) = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window,[], feedf.method, 0, 0);
                
                
                % ibi0
                % -------------------------------------------
                % STORE IN LONG & WIDE FORMAT
                % -------------------------------------------
                
                
                
                % -------------------------------------------
                % STORE IN LONG FORMAT
                % -------------------------------------------
                
                % FACTORS
                longF{i,1} = suji ; %subject id
                longF{i,2} = roi_i; % roi
                longF{i,3} = roilabels{roi_i}; % roi
                longF{i,4} = condi; % condition (1 = low; 2 = hi)
                % MEASURES
                longF{i,5} = round(outputP.peakminusbaseline(roi_i, condi),2); % outputP(roi_i, condi)
                longF{i,6} = round(outputP.wmean(roi_i, condi),2);%
                longF{i,7} = round(outputP.onsetnoise_ms(roi_i, condi),2);%
                longF{i,8} = output_abs(roi_i, condi).risingthreshlat; %
                
                %                     longF(i,6) = round(NEWMEASURE);%
                %                     longF(i,5) = outputP.peaklatency; %
                
                
                % -------------------------------------------
                % STORE IN WIDE FORMAT (for grouped barplots)
                % -------------------------------------------
                wideF.peakminusbasel(roi_i, condi, suji) = longF{i,5};
                wideF.wmean(roi_i, condi, suji) =  longF{i,6};%
                wideF.onsetnoise_ms(roi_i, condi, suji) =  longF{i,7};%
                wideF.risingthreshlat(roi_i, condi, suji) = longF{i,8}; %
                
                
                i = i+1;
                
            end % roi_i
            clear  roiwave movieave  sel_trials
            
        end % condi
        
        % SUB-LOOP 2B:  NORMALIZATION
        
        
        for condi = 1:n
            
            for roi_i = makeRow(selroi)
                
                roiwave =roiwaves(:,roi_i, condi, suji); % GET BACK FROM PREVIOUS LOOP
                
                % STORE NORMALIZED MEASURES
                longF{j,9} =round(longF{j,5}/max(outputP.peakminusbaseline(:)),2); %
                longF{j,10} = round(longF{j,6}/max(outputP.wmean(:)),2); %
                
                % -------------------------------------------
                % STORE IN WIDE FORMAT (for grouped barplots)
                % -------------------------------------------
                % .. . AND NORMALIZE AS WELL
                wideF_norm.peakminusbasel(roi_i, condi, suji) = longF{j,9};
                wideF_norm.wmean(roi_i, condi, suji) =  longF{j,10};%
                
                j = j+1;
                
            end % roi_i
            
            
        end % condi
        
        clear  output_abs outputP
        
        
    end % suji
    
    
    %% THIRD 'suji' LOOP: PLOT AND BAR-PLOTS
    
    for suji  =  1:size(groupplot,1)
        nfish = groupplot{suji,1};
        VSDI = TORus('load', nfish); %we need it to get the info of the fish for plot title
        
        %----------------------------------------------------------------
        % SELECT ROI
        %----------------------------------------------------------------
        
        switch roikind
            case 'circle'
                selroi =name2idx(selroinames, VSDI.roi.labels_circ);
                roilabels = VSDI.roi.labels_circ;
            case 'anat'
                selroi =name2idx(selroinames, VSDI.roi.labels);
                roilabels = VSDI.roi.labels;
        end
        
        %     for condi = 1:n
        %
        %         for roi_i = makeRow(selroi)
        %             %
        %             roiwave =roiwaves(:,roi_i, condi, suji); % GET BACK FROM PREVIOUS LOOP
        
        %% PLOT FOR EACH SUBJECT : BEFORE CLOSING THE THIRD 'for suji' LOOP
        % --------------------------------------------------------
        % PLOTS : HARDWIRED FOR 2-ROIS AND 2-CONDITIONS !!! plot both rois from both conditions (in each
        % subject)
        % --------------------------------------------------------
        if ploton
            figure
            ax(1) = subplot(2,3,[1 2 4 5]);
            plot(VSDI.timebase(idxrange), roiwaves(:,selroi(1), 1, suji), '--b', 'linewidth', 1, 'DisplayName',strcat(selroinames{1}, '-low'))
            hold on
            plot(VSDI.timebase(idxrange), roiwaves(:,selroi(1), 2, suji), '-b', 'linewidth', 1, 'DisplayName',strcat(selroinames{1}, '-high'))
            
            plot(VSDI.timebase(idxrange), roiwaves(:,selroi(2), 1, suji), '--r',  'linewidth', 1, 'DisplayName',strcat(selroinames{2}, '-low'))
            plot(VSDI.timebase(idxrange), roiwaves(:,selroi(2), 2, suji), '-r',  'linewidth', 1, 'DisplayName',strcat(selroinames{2}, '-high'))
            
            if flagblank
                plot(VSDI.timebase(idxrange), roiwaves(:,selroi(1), 3, suji), ':b', 'linewidth', 1, 'DisplayName',strcat(selroinames{1}, '-control'))
                plot(VSDI.timebase(idxrange), roiwaves(:,selroi(2), 3, suji), ':r', 'linewidth', 1, 'DisplayName',strcat(selroinames{2}, '-control'))
            end
            
            xlim([VSDI.timebase(idxrange(1)), VSDI.timebase(idxrange(end))])
            legend()
            tname = [num2str(VSDI.ref) '(' ref_movie(2:end) '-' roikind ')' ];
            
            ax(2) = subplot(2,3,3);
            
            switch roikind
                case 'circle'
                    roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R, ax(2));
                case 'anat'
                    roi_preview_multiple(VSDI.crop.preview,  VSDI.roi.manual_poly(selroi,:), ax(2))
            end
            ax(2).Visible = 'off';
            %     title([num2str(VSDI.ref) 'roi preview:' selroinames{:}])
            
            
            sgtitle(tname)
            saveas(gcf, fullfile(savein, ['plot' tname  '_reject' num2str(reject_on) refcase '.jpg']), 'jpg')
            close
            %                 --------------------------------------------------------
        end %ploton
        
        
        %         end %for roi_i
        %     end % for condi
        
        % --------------------------------------------------------
        % BARPLOTS : HARDWIRED FOR 2-ROIS AND 2-CONDITIONS !!! plot both rois from both conditions (in each
        % subject)
        % --------------------------------------------------------
        if plotbar_on
            
            sgtitle(['fish-' num2str(VSDI.ref)])
            subplot(2,3,1)
            bar(wideF.peakminusbasel(selroi,:, suji))
            set(gca,'xticklabel', selroinames )
            
            title('peak')
            
            subplot(2,3,2)
            bar(wideF.wmean(selroi, :, suji)); %
            set(gca,'xticklabel', selroinames )
            title('wmean')
            
            subplot(2,3,3)
            bar(wideF.onsetnoise_ms(selroi, :, suji)); %
            title('onset noise ms')
            set(gca,'xticklabel', selroinames )
            
            subplot(2,3,4)
            bar(wideF.risingthreshlat(selroi, :, suji)); %
            title('onset latency 20%(abs)')
            set(gca,'xticklabel', selroinames )
            
            
            %                    subplot(2,3,5)
            %                    bar(wideF.risingthreshlat(selroi, :)); %
            %                    title('JUST FOR LEGEND')
            %                    legend( {'lo','hi','con'})
            
            
            subplot(2,3,5)
            bar(wideF_norm.peakminusbasel(selroi,:, suji))
            set(gca,'xticklabel', selroinames )
            title('peak (norm)')
            
            subplot(2,3,6)
            bar(wideF_norm.wmean(selroi, :, suji)); %
            set(gca,'xticklabel', selroinames )
            title('wmean(norm)')
            
            saveas(gcf, fullfile(savein, ['barplots' num2str(VSDI.ref) 'reject' num2str(reject_on) refcase '.jpg']), 'jpg')
            close
        end
    end % for suji
    
    %% -------------------------------------------
    % EXPORT FOR R
    % -------------------------------------------
    excelname = fullfile(savein, ['long_format_forR_group.xls']);
    labels = {'id' 'roi' 'roi n' 'cond' 'peak' 'wmean' 'noise onset' 'onset 20abs_ms' 'peak(norm)' 'wmean(norm)'};
    
    for col = 1:numel(labels)
        longF{1,col}= labels{col};
    end
    
    
    % write output (new sheet for each fish
    writecell (longF, excelname, 'sheet',  [roikind ref_movie 'rej' num2str(reject_on) refcase])
    
    writecell (labels, excelname, 'sheet', 'labels')
    writecell (groupplot, excelname, 'sheet','group')
    writecell (params, excelname, 'sheet','param')
    
    blob()
end % for reject_on
