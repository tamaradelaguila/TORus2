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
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new.mat') % ATTENTION: select cases manually (there are repetitions)

% ///////////////////////////////////////////////////////////
% SETTINGS

% selroinames = {'dm4m_R',  'dm2_R'};
% selroinames = {'dm4m_R2',  'dm2_R2'};

selroinames = {'dm4m_R',  'dm2_R', 'dm3_R' ,'dm1_R','dldm_R'};%dm3,
% selroinames = {'dm4m_R2',  'dm2_R2', 'dm3_R2' ,'dm1_R','dldm_R2'};%dm3,
% selroinames = {'dm4m_R',  'dm4_L', 'dm2_R' ,'dm2_L'};
% selroinames = {'dm4m_R2','dm2_R2','dm1_R','dldm_R2'}; %4roi
% selroinames = {'dm4m_R2','dm2_R2','dm1_R','dldm_R2','dm4m_L2','dm2_L2','dm1_L','dldm_L2'}; %8roi
% selroinames = {'dm4m_R2','dm2_R2','dm4m_L2','dm2_L2'}; %4roi

roikind = 'circle'; %
% roikind = 'anat';

dataunits = 'dF'; % '%F' 'dF'
switch dataunits
    case 'dF'
        ref_movie = '_18filt6'; % input movie '_17filt5'
    case '%F'
        ref_movie = '_21filt6'; 
end

analysisref = 'new_group6_RECOV'; % MANUALLY SET!!! extra info for the name. Set group according to the rows of 'groupplot' selected in;   for suji  =  [1 3 4 9] lateral_group7_RECOV
reject_on = 3;

onset_factor = 0.2; % of max value to set rising threshold (to onset latency)

plotmaps = 0;
savemaps = 0; % if plotmaps =1
getR = 1;
get_excel_indiv = 0; %not working at the moment

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/Zscore3_percF' ;%@ SET

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

if strcmpi( analysisref, 'new_group2_RECOV')  || strcmpi( analysisref, 'lateral_group2_RECOV') 
    sel_subjects = [1:4,9];
elseif strcmpi( analysisref, 'new_group3_RECOV') || strcmpi( analysisref, 'lateral_group3_RECOV') 
    sel_subjects = [5:8];
elseif strcmpi( analysisref, 'new_group4_RECOV') || strcmpi( analysisref, 'lateral_group4_RECOV') 
    sel_subjects = [1 3 4 9];
elseif strcmpi( analysisref, 'new_group6_RECOV') || strcmpi( analysisref, 'lateral_group6_RECOV') 
    sel_subjects = [2 3 8 9 11];
elseif strcmpi( analysisref, 'new_group7_RECOV') || strcmpi( analysisref, 'lateral_group7_RECOV')
    sel_subjects = [2 8 9 11];
end

    %----------------------------------------------------------------
    % @SET: REJECT SETTINGS
    %----------------------------------------------------------------
    
    % Subsettings:
    setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET+
    setting.forcein = 0; %@ SET
    

    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    i = 2; % counter for long format rows: the first will be the labels
    si = 0; % counter for subjects list
    
    for suji  =  sel_subjects %1:size(groupplot,1) SELECT included fish+condition
        
        si = si+1;
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
rejectidx  = compute_rejectidx(VSDI, reject_on, setting);
        
        
        
%         params{8,1} = 'rejected';

%         params{8,1} = rejectidx';

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
        % FIRST LOOP THROUGH CONDITIONS AND PIXELS TO GET PIXEL-MEASURES
        % -------------------------------------------
        
        ci = 0;
        for condition =  makeRow(cond_lohi)
            ci = ci+1;
            
            cond_blank = force0ending(condition);
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
            % CALCULATE MEASURES for each pixel and condition (and
            % store in 'maps')
            % -------------------------------------------
                    
                    for xi = 1:size(movieave,1)
                        for yi = 1:size(movieave,2)
                            pixelwave = movieave(xi, yi,:);
                            
                            temp = devo_peak2peak(pixelwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
                            
                            maps.peak(xi,yi,ci) = temp.peakminusbasel;
                            maps.wmean(xi,yi,ci) = temp.wmean;
%                             maps.peaklat(xi,yi,ci) = temp.peaklatency;
%                             % pixel-wise peak latency
                            
                            %slope
                                idx0 = slope.windowidx(1);
                                idxend = slope.windowidx(end);
                                waveW = pixelwave(idx0:idxend);
                                slopemean = mean(diff(waveW));
                                maps.slope(xi,yi,ci) = slopemean;
                                clear pixelwave waveW slopemean

                        end % for yi
                    end % for xi
        
            % get maxval for ONSET LATENCY THRESHOLD
             maps.max(ci) = max(movieave(:));

             
             %---------------------------------------------------------------
             ... MEASURES FROM WAVES:  SECOND LOOP THROUGH CONDITIONS AND PIXELS TO GET ONSET LATENCY
             ...(once the threshold has been set respect to the maximum values
             % ---------------------------------------------------------------
% ploti = 0;
for roi_i = makeRow(selroi)
%     ploti = ploti+1;
%     subplot(1,length(selroi), ploti)
    roimask= masks(:,:,roi_i);
    
    roiwave =  roi_TSave(movieave,roimask);
    
    temp = devo_peak2peak(roiwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
    peaklat(roi_i,ci) = temp.peaklatency;
    waves.peak(roi_i,ci) = temp.peakminusbasel;
    waves.mean(roi_i,ci) = temp.wmean;
    
    %slope
    idx0 = slope.windowidx(1);
    idxend = slope.windowidx(end);
    waveW = roiwave(idx0:idxend);
    slopemean = mean(diff(waveW));
    waves.slope(roi_i,ci)  = slopemean;
    clear pixelwave waveW slopemean
    
% %     % PLOT TO CHECK
%     plot(timebase_adj, roiwave, 'linewidth', 1.3)
%     xline(temp.peaklatency)
%     title([num2str(VSDI.ref) ':' roilabels{roi_i}])
%     pause
    
    clear temp roiwave
    
end
% sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Peaklatency from wave .Cond:' num2str(condition)])
% 
% localname = ['Peaklat_from_avewave_cond' num2str(condition)] ;
% savename = [num2str(VSDI.ref) '_' ref_movie '_rej' num2str(reject_on) '_' ];
% 
% saveas(gcf, fullfile(savein, [savename localname '.jpg']), 'jpg') % prints it as you see them
% close

        end %for condition
        
        % ---------------------------------------------------------------
        % LATENCY FROM MAPS:  SECOND LOOP THROUGH CONDITIONS AND PIXELS TO GET ONSET LATENCY
        % (once the threshold has been set respect to the maximum values
        % ---------------------------------------------------------------

           feedf.window.risingthresh = onset_factor*max(maps.max);
 
           ci = 0;
           for condition =  makeRow(cond_lohi)
               
               ci = ci+1;
               
               % ----------------------------------------------------
               % SELECT CASES  AND AVERAGE MOVIE
               % ----------------------------------------------------
               sel_trials = find(VSDI.condition(:,1)==condition);
               
               if reject_on
                   sel_trials = setdiff(sel_trials, rejectidx);
               end
               
               % --------------------------------------------------------------------------
               % GET AVERAGE MOVIE. Use timerange set
               % --------------------------------------------------------------------------
               
               movieave = mean(movies(:,:,idxrange,sel_trials),4);
               movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
               
               for xi = 1:size(movieave,1)
                   for yi = 1:size(movieave,2)
                       pixelwave = squeeze(movieave(xi, yi,:));
                       
                       temp = devo_peak2peak(pixelwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
                       maps.onsetlat(xi,yi,ci) = temp.risingthreshlat;
                       clear pixelwave temp
                   end
               end
               
           end


        %% SPATIAL Z-SCORE AMONG CONDITIONS (for each measure)
        
%         dim = size(maps.peak);
%         temp_peak = reshape(maps.peak, [dim(1)*dim(2) dim(3)]);
%         temp_wmean = reshape(maps.wmean, [dim(1)*dim(2) dim(3)]);
%         temp_slope = reshape(maps.slope, [dim(1)*dim(2) dim(3)]);
%         temp_onsetlat = reshape(maps.slope, [dim(1)*dim(2) dim(3)]);
        
        Zpeak = zscore(maps.peak, 0, 'all');
        Zwmean = zscore(maps.wmean, 0, 'all');
        Zslope = zscore(maps.slope, 0, 'all');
        
        mapsZ.peak = Zpeak;
        mapsZ.wmean = Zwmean;
        mapsZ.slope = Zslope;

        clear Zpeak Zmean Zslope
%         mapsZ.peak =reshape(Zpeak, [dim(1) dim(2) dim(3)]);
%         mapsZ.wmean =reshape(Zwmean, [dim(1) dim(2) dim(3)]);
%         mapsZ.slope = reshape(Zslope, [dim(1) dim(2) dim(3)]);

        
    
%%         %% -------------------------------------------
%         % GET INDIVIDUAL  MEASURES OF ROI WAVES FOR EACH TRIAL - TO DEBUG
%         % -------------------------------------------
%         
%         if get_excel_indiv
%             ci = 0;
%             ti = 1; ri = 3;
%             for condition =  makeRow(cond_lohi)
%                 ci = ci +1;
%                 
%                 sel_trials = find(VSDI.condition(:,1)==condition);
%                 if reject_on
%                     sel_trials = setdiff(sel_trials, rejectidx);
%                 end
%                 
%                 for triali = makeRow(sel_trials) %each trial will be saved in a row
%                     ti = ti+1;
%                     movie2plot = movies(:,:,idxrange,triali);
%                     ri= 3;
%                     for roii = makeRow(selroi) %each roi in a column
%                         ri= ri+1;
%                         roimask = masks(:,:, roii);
%                         %                             meanF0 = squeeze(mean(VSDmov.F0(:,:,triali),3));
%                         roiwave = roi_TSave(movie2plot,roimask);
%                         temp = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window, [], feedf.method, 0, 0);
%                         
%                         % trial identification / kind
%                         trialmeasure.peak{ti,1} = triali;
%                         trialmeasure.peak{ti,2} = VSDI.trialref(triali);
%                         trialmeasure.peak{ti,3} = condition;
%                         % measures: PEAK
%                         trialmeasure.peak{ti,ri} = round(temp.peakminusbasel,2);
%                         
%                         % trial identification / kind
%                         trialmeasure.wmean{ti,1} = triali;
%                         trialmeasure.wmean{ti,2} = VSDI.trialref(triali);
%                         trialmeasure.wmean{ti,3} = condition;
%                         % measures: WMEAN
%                         trialmeasure.wmean{ti,ri} = round(temp.wmean,2);
%                         
%                         clear roimask roiwave meanF0
%                         labels{1,1} = 'idx';
%                         labels{1,2} = 'trial';
%                         labels{1,3} = 'condition';
%                         labels{1,ri} = roilabels{roii};
%                         
%                         % TRIALS TABLE TO LATER REJECT BASED ON STD
%                         % trial identification / kind
%                         trialmat.peak(ti,1) = triali;
%                         trialmat.peak(ti,2) = VSDI.trialref(triali);
%                         trialmat.peak(ti,3) = condition;
%                         % measures: PEAK
%                         trialmat.peak(ti,ri) = round(temp.peakminusbasel,2);
%                         
%                         % trial identification / kind
%                         trialmat.wmean(ti,1) = triali;
%                         trialmat.wmean(ti,2) = VSDI.trialref(triali);
%                         trialmat.wmean(ti,3) = condition;
%                         % measures: WMEAN
%                         trialmat.wmean(ti,ri) = round(temp.wmean,2);
%                         
%                     end % forroii
%                 end % for triali
%             end %for condition
%             % add labels
%             trialmeasure.peak(1,1:length(labels)) = labels(1:end);
%             trialmeasure.wmean(1,1:length(labels)) = labels(1:end);
%             
%             excelname = [analysisref 'Zscore_individual_trials' ref_movie '_rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls'];
%             labels = {'idx' 'trial' 'cond'};
%             
%             % write output (new sheet for each fish
%             writecell (trialmeasure.peak, fullfile(savein , excelname), 'sheet',  [num2str(VSDI.ref) 'peak'])
%             writecell (trialmeasure.wmean, fullfile(savein,  excelname), 'sheet',  [num2str(VSDI.ref) 'wmean'])
%             
%             clear trialmeasure
%             
%             
%         end %if get_excel_indiv
        

        %% EXTRACT ROI MEASURES AND STORE IN LONG FORMAT (for R)
        
        ci = 0;
        for condition =  makeRow(cond_lohi)
            ci = ci +1;
            
            if getR
                % -------------------------------------------------------
                % CALCULATE MEASURE FOR EACH ROI AND STORE IN LONG FORMAT
                % -------------------------------------------------------
                for roi_i = makeRow(selroi)
                    
                    roimask= masks(:,:,roi_i);
                    roipeak = sum(maps.peak(:,:,ci).*roimask) / sum(roimask) ;
                    roiwmean = sum(maps.wmean(:,:,ci).*roimask) / sum(roimask) ;
                    roislope = sum(maps.slope(:,:,ci).*roimask) / sum(roimask) ;
                    
                    roipeakZ = sum(mapsZ.peak(:,:,ci).*roimask) / sum(roimask) ;
                    roiwmeanZ = sum(mapsZ.wmean(:,:,ci).*roimask) / sum(roimask) ;
                    roislopeZ =  sum(mapsZ.slope(:,:,ci).*roimask) / sum(roimask) ;
                    
                    roionsetlat = sum(maps.onsetlat(:,:,ci).*roimask, 'omitnan') / sum(roimask) ;
                    roipeaklat =  peaklat(roi_i,ci);

                    % FACTORS
                    longF{i,1} = suji ; %subject id
                    longF{i,2} = roi_i; % roi
                    longF{i,3} = roilabels{roi_i}; % roi
                    longF{i,4} = ci; % condition (1 = low; 2 = hi; 3 = blank)
                    % MEASURES
                    longF{i,5} = round(roipeak,2); % outputP(roi_i, condi)
                    longF{i,6} = round(roiwmean,2);%
                    longF{i,7} = round(roislope,2);%
                    longF{i,8} = round(roionsetlat);%

                    longF{i,9} = round(roipeakZ,2); % outputP(roi_i, condi)
                    longF{i,10} = round(roiwmeanZ,2);%
                    longF{i,11} = round(roislopeZ,2);%
                    longF{i,12} = round(roipeaklat);


                    % set hemisphere factor
                    if strcmpi( roilabels{roi_i}(end), 'R') || strcmpi(roilabels{roi_i}(end-1),'R')
                        longF{i,13} = 'R';%
                    elseif strcmpi( roilabels{roi_i}(end), 'L') || strcmpi(roilabels{roi_i}(end-1),'L')
                        longF{i,13} = 'L';%
                    end
                    
                    longF{i,14} = roilabels{roi_i}(1:4);%in case we want to perform anova roi*hemisf
                    
                    longF{i,15} = round(waves.peak(roi_i,ci),3); % outputP(roi_i, condi)
                    longF{i,16} = round(waves.mean(roi_i,ci),3);%
                    longF{i,17} = round(waves.slope(roi_i,ci),5);%


                    i = i+1;
                    clear roimask  roipeak  roiwmean roislope roipeakZ roiwmeanZ roislopeZ
                end %for roi_i
            end % for if getR
        end % for condition
        
        %% PLOT MAPS
        if plotmaps
            
            savename = [num2str(VSDI.ref) '_' ref_movie '_rej' num2str(reject_on) '_'];
            %
            ncond = length(cond_lohi);
            
            %             % ------------------------------------------------------------------
            %             % PLOT MAPS OF AVERAGE MEASURES CONDITION-WISE
            %             % ------------------------------------------------------------------
            %             figure
            %             ci = 0;
            %             for condition =  makeRow(cond_lohi)
            %                 % peak in the first row
            %                 ci = ci+1;
            %                 ax(ci) = subplot(2,ncond,ci);
            %                 imagesc(maps.peak(:,:,ci))
            %                 axis image
            %                 colorbar
            %                 %
            %                 set(gca, 'clim', [0 max(maps.peak(:))*0.8])
            %                 condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
            %                 tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            %                 title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
            %
            %                 % wmean in the second row
            %                 ax(ci+ncond) = subplot(2,ncond,ci+ncond);
            %                 imagesc(maps.wmean(:,:,ci))
            %                 axis image
            %                 colorbar
            %                 set(gca, 'clim', [0 max(maps.wmean(:))*0.8])
            %
            %
            %             end
            %
            %             sgtitle([num2str(VSDI.ref), '(',ref_movie,  ')', 'up: peak; down: onsetA'])
            %             localname = ['MAPS'] ;
            %
            %             if savemaps
            %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
            %                 close
            %             end
            
%             % ------------------------------------------------------------------
%             ... PLOT MAPS OF SPATIAL Z-SCORE OF AVERAGE MEASURES CONDITION-WISE
%             % ------------------------------------------------------------------
            
            % colormap like jet but changing initial values 
            ccmap = jet;
            % remove darker colors: 
            darkblue = ccmap(4,:);
            ccmap = removerows(ccmap,'ind',[1:8, 1:3]);

            % change first color and interpolate
            ccmap(1,:) = darkblue;
            flag = 3;
            R = linspace(ccmap(1,1), ccmap(flag,1), flag);
            G = linspace(ccmap(1,2), ccmap(flag,2), flag);
            B = linspace(ccmap(1,3), ccmap(flag,3), flag);
            ccmap(1:flag,:) = [R; G; B]'; 
            
            cclim = [0 4];
%             %         ccmap(1,:) = [0 0 0];
%             
% 
%             % PLOT PEAK + WMEAN (in the same figure)
%             % ------------------------------------------------------------------
%             figure
%             ci = 0;
%             for condition =  makeRow(cond_lohi)
%                 % peak in the first row
%                 ci = ci+1;
%                 ax(ci) = subplot(2,ncond,ci);
%                 
%                 im1 = mapsZ.peak(:,:,ci);
%                 %                 im1 = interp2(im1, 5, 'nearest');
%                 %                 im1(~VSDI.crop.preview) = 0;
%                 
%                 alphamask = ones(size(im1))*0.6;
%                 
%                 %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
%                 %                 alphamask =  interp2(alphamask, 5, 'nearest');
%                 
%                 back = VSDI.backgr(:,:,VSDI.nonanidx(1));
%                 back = interp2(back,5);
%                 %                 imagesc(im1)
%                 %                 axis image
%                 %                 set(gca, 'clim',cclim)
%                 %                 colorbar; colormap(ccmap)
%                 %
%                 plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
%                 ax(ci).Visible = 'off';
%                 % STOPS THE CODE FOR CHECKING ROI CENTERS
%                 %             if ci ==2
%                 %                 return
%                 %             end
%                 
%                 condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
%                 tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
%                 title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
%                 
%                 % wmean in the second row
%                 ax(ci+ncond) = subplot(2,ncond,ci+ncond);
%                 
%                 im2 = mapsZ.wmean(:,:,ci);
%                 %                 im2(~VSDI.crop.preview) = 0;
%                 
%                 alphamask = ones(size(im2))*0.6;
%                 %                 alphamask= im2  >0;
%                 
%                 back = VSDI.backgr(:,:,VSDI.nonanidx(1));
%                 back = interp2(back,5);
%                 plot_framesoverlaid2(im2, back, alphamask, 0, ax(ci+ncond), cclim, 1 , 0, ccmap)
%                 ax(ci+ncond).Visible = 'off';
%                 
%             end
%             
%             sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Zmap. up: peak; down: onsetA. Cond:' num2str(cond_lohi)])
%             localname = ['MAPS_Zscore'] ;
%             
%             if savemaps
%                 
%                 %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%                 print(fullfile(savein, [analysisref savename localname '.jpg']),'-r900','-djpeg') % prints it as you see them
%                 
%                 close
%             end
%             
%             
%             % PLOT PEAK 
%             % ------------------------------------------------------------------
%             figure
%             ci = 0;
%             for condition =  makeRow(cond_lohi)
%                 % peak in the first row
%                 ci = ci+1;
%                 ax(ci) = subplot(1,ncond,ci);
%                 
%                 im1 = mapsZ.peak(:,:,ci);
%                 %                 im1 = interp2(im1, 5, 'nearest');
%                 %                 im1(~VSDI.crop.preview) = 0;
%                 
%                 alphamask = ones(size(im1))*0.6;
%                 
%                 %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
%                 %                 alphamask =  interp2(alphamask, 5, 'nearest');
%                 
%                 back = VSDI.backgr(:,:,VSDI.nonanidx(1));
%                 back = interp2(back,5);
%                 %                 imagesc(im1)
%                 %                 axis image
%                 %                 set(gca, 'clim',cclim)
%                 %                 colorbar; colormap(ccmap)
%                 %
%                 plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
%                 ax(ci).Visible = 'off';
%                 % STOPS THE CODE FOR CHECKING ROI CENTERS
%                 %             if ci ==2
%                 %                 return
%                 %             end
%                 
%                 condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
%                 tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
%                 title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
%                 
%                 
%             end
%             
%             sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'PEAK-Z map. Cond:' num2str(cond_lohi)])
%             localname = ['MAPS_peakZ'] ;
%             
%             if savemaps
%                 
%                 %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%                 print(fullfile(savein, [analysisref savename localname '.jpg']),'-r900','-djpeg') % prints it as you see them
%                 
%                 close
%             end% 

            % PLOT WMEAN
            % ------------------------------------------------------------------
            figure
            ci = 0;
            for condition =  makeRow(cond_lohi)
                % peak in the first row
                ci = ci+1;
                
                condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
                tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
                title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
                
                % wmean in the second row
                ax(ci) = subplot(1,ncond,ci);
                
                im2 = mapsZ.wmean(:,:,ci);
                %                 im2(~VSDI.crop.preview) = 0;
                
                alphamask = ones(size(im2))*0.6;
                %                 alphamask= im2  >0;
                
                back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                back = interp2(back,5);
                plot_framesoverlaid2(im2, back, alphamask, 0, ax(ci), cclim, 1 , 0, ccmap)
                ax(ci).Visible = 'off';
                
            end
            
            sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'onsetA-Z. Cond:' num2str(cond_lohi)])
            localname = ['MAPS_wmeanZ'] ;
            
            if savemaps
                %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
                
                close
            end
            
%             % PLOT SLOPE
%             % ------------------------------------------------------------------
%             figure
%             ci = 0;
%             for condition =  makeRow(cond_lohi)
%                 % peak in the first row
%                 ci = ci+1;
%                 ax(ci) = subplot(1,ncond,ci);
%                 
%                 im1 = mapsZ.slope(:,:,ci);
%                 %                 im1 = interp2(im1, 5, 'nearest');
%                 %                 im1(~VSDI.crop.preview) = 0;
%                 
%                 alphamask = ones(size(im1))*0.6;
%                 
%                 %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
%                 %                 alphamask =  interp2(alphamask, 5, 'nearest');
%                 
%                 back = VSDI.backgr(:,:,VSDI.nonanidx(1));
%                 back = interp2(back,5);
%                 %                 imagesc(im1)
%                 %                 axis image
%                 %                 set(gca, 'clim',cclim)
%                 %                 colorbar; colormap(ccmap)
%                 %
%                 plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
%                 ax(ci).Visible = 'off';
%                 % STOPS THE CODE FOR CHECKING ROI CENTERS
%                 %             if ci ==2
%                 %                 return
%                 %             end
%                 
%             end
%             
%             sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Zmap. Slope. Cond:' num2str(cond_lohi)])
%             localname = ['Zscore_slope'] ;
%             
%             if savemaps
%                 
%                 %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%                 print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
%                 
%                 close
%             end
            
%             % PLOT ONSET LATENCY
%             % ------------------------------------------------------------------
%             figure
%             ci = 0;
%             for condition =  makeRow(cond_lohi)
%                 % peak in the first row
%                 ci = ci+1;
%                 ax(ci) = subplot(1,ncond,ci);
%                 
%                 im1 = maps.onsetlat(:,:,ci);
%                 %                 im1 = interp2(im1, 5, 'nearest');
%                 %                 im1(~VSDI.crop.preview) = 0;
%                 
%                 alphamask = ones(size(im1))*0.6;
%                 
%                 %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
%                 %                 alphamask =  interp2(alphamask, 5, 'nearest');
%                 
%                 back = VSDI.backgr(:,:,VSDI.nonanidx(1));
%                 back = interp2(back,5);
%                 %                 imagesc(im1)
%                 %                 axis image
%                 %                 set(gca, 'clim',cclim)
%                 %                 colorbar; colormap(ccmap)
%                 %
%                 
%                 flipmap =flipud(ccmap);
%                 plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), [50 80], 1, 0, flipmap)
%                 ax(ci).Visible = 'off';
% 
%                 % STOPS THE CODE FOR CHECKING ROI CENTERS
%                 %             if ci ==2
%                 %                 return
%                 %             end
%             end
%             
%             sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Onsetlat' num2str(onset_factor*100) '%. Cond:' num2str(cond_lohi)])
%             localname = ['onsetlat' num2str(onset_factor*100) '_'] ;
%             
%             if savemaps
%                 
%                 % saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%                 print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
%                 
%                 close
%             end
%             
% 
            % ------------------------------------------------------------------
            % PLOT ROIS USED
            % ------------------------------------------------------------------
            figure
            ax1 = subplot(1,1,1); 
            switch roikind
                case 'circle'
                    centers = VSDI.roi.circle.center(selroi, :) ;
                    roicirc_preview_multiple(VSDI.crop.preview, centers, VSDI.roi.circle.R, ax1);
                    
                case 'anat'
                    roi_preview_multiple(VSDI.crop.preview, VSDI.roi.manual_poly(selroi,:), ax1);
            end
            
            localname = [num2str(VSDI.ref), ref_movie,'roipreview',num2str(numel(selroi)) ];
            
            if savemaps
                saveas(gcf, fullfile(savein, [analysisref localname '.jpg']), 'jpg')
                close
            end
        end % if plotmaps
        
        
        %         % PLOT PEAK VALUES
        %         ax2 = subplot(1,3,2);
        %         temp = barplot.peak(suji,:,selroi)';
        %         legend(selroinames)
        %         title(num2str(suji))
        %
        %         localref = [num2str(VSDI.ref) '- cond=' num2str(condition) '-' ref_movie '-' refcase  '-reject' num2str(reject_on) ];
        %         sgtitle(localref)
        %
        %         savename = [ localref   '-' num2str(numel(selroi)) roikind 'ROI - PART2 WAVES'];
        % %         localname = 'tilewaves'];
        %         saveas(gcf, fullfile(savein, [savename localname '.jpg']), 'jpg')
        %         close
        
        %
        clear maps mapsZ
        
        groupplot_print(suji,:) = groupplot(suji,1:end);
    end % for suji
    
    %% -------------------------------------------
    % EXPORT FOR R
    % -------------------------------------------
    if getR
        excelname = fullfile(savein, [analysisref '_' dataunits '_Zscore_long_format_forR_group' ref_movie '_rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls']);
        labels = {'id' 'roi' 'roi n' 'cond' 'peak' 'wmean' 'slope' 'onsetlat' 'Zpeak' 'Zmean' 'Zslope' 'peaklat' 'hemis', 'roii',  'wavepeak' 'wavemean' 'waveslope' };
        
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
% 27/02/22 - Add 'waves.peak', 'waves.mean', 'waves.slope': measures 
... straight from roiwave (extracted from the average movie), instead of from Zmaps. To test
... whether %F needs some extra normalization or not.
% 16/02/22 - Add onset-latency measure and delete 'switch refcase' (so only
% the normal dF condition is considered). To use blank-substraction
% methods, go to the older versions of the code
