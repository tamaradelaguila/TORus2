%% GET WAVES 

% z-spatial all  conditions for each fish
% extract roi from z-spatially maps >>> definite measures
close all
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus';
user_settings
cd(W)

%----------------------------------------------------------------
% @SET: fish +  conditions
%---------------------------------------------------------------
nfish = 14;
 condition = [404];

% movie to preview wave
ref_movie= '_18filt6';% '_15filt5', '_17filt5' ; '_18filt6'

% selroinames = {'dm4m_R',  'dm2_R'};
% selroinames = {'dm4m_R2',  'dm2_R2'};

selroinames = {'dm4m_R',  'dm2_R', 'dm3_R' ,'dm1_R','dldm_R'};%dm3,
nHemis = 1;

roikind = 'circle'; %
% roikind = 'anat';

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis

% CONFIGURATION PARAMETERS
% -------------------------------------------
measure_mehod =  '%F_meth3_avewise_roi';

analysisref = 'new4_group10_pulsito'; % MANUALLY SET!!! extra info for the name. Set group according to the rows of 'groupplot' selected in;   for suji  =  [1 3 4 9] lateral_group7_RECOV
reject_on = 0; %before 10/09/22, reject_on = 3;

savewaves = 0;

%% --- 

switch measure_mehod
    case 'dF'
        ref_movie = '_18filt6'; % input movie '_17filt5'
dataunits = 'dF'; % '%F' 'dF'
    case '%F_meth1_from%movie' % from movie
        ref_movie = '_21filt6';
        dataunits = '%F'; % '%F' 'dF'

    case  '%F_meth2_trialwise'
        ref_movie = '_18filt6'; % input movie '_17filt5'
        dataunits = '%F'; % '%F' 'dF'
        
    case '%F_meth3_avewise1' % first
        ref_movie = '_18filt6'; % input movie '_17filt5'
        dataunits = '%F'; % '%F' 'dF'

end


        % -------------------------------------------------------
        % SET COLORMAP
        % -------------------------------------------------------
            cmap = roicolors();
        if nHemis ==1
            cmap = cmap(1:2:end,:); %roicolors map has double values for 2 hemispheres
        end


%% COMMON TO ALL METHODS

%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------

VSDI = TORus('load', nfish);
VSDmov = TORus('loadmovie',nfish,ref_movie);
movies = VSDmov.data ;

    %----------------------------------------------------------------
% GET INDEXES OF TIMERANGE
%----------------------------------------------------------------
idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values

timebase_adj = VSDI.timebase(idxrange);

% ----------------------------------------------------------------
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


%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

% Subsettings:
setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET+
setting.forcein = 0; %@ SET

% ----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------

rejectidx = [];
rejectidx  = compute_rejectidx(VSDI, reject_on, setting);

%%  -------------------------------------------
% METHOD 2
% -------------------------------------------
switch measure_mehod
    case '%F_method2'

    %         cond_blank = force0ending( condition); %UNUSED when there is a blank trial for each  condition

    % -------------------------------------------
    % SELECT CASES  
    % -------------------------------------------
    sel_trials = find(VSDI. condition(:,1)== condition);

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
% CALCULATE %F WAVES FOR EACH ROI, AND ITS CORRESPONDING MEASURE
% -------------------------------------------

% GET %F AVE FOR EACH TRIAL AND ROI AND STORE FOR LATER
% AVERAGING
t = 0;
for ti = makeRow(sel_trials)
    t = t+1;
    movie = movies(:,:,idxrange,ti);
    movie(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave

    F0 = VSDmov.F0(:,:,ti);

    for roi_i = makeRow(selroi)
        %     ploti = ploti+1;
        %     subplot(1,length(selroi), ploti)
        roimask= masks(:,:,roi_i);
        roiwave =  roi_TSave_percF_roiwise(movie,roimask, F0);
        percFwave(:,roi_i,t) = roiwave;
        clear roimask roiwave

    end % for roi

    clear movie
end % ti (seltrials)


% GET MEASURES ROI-WISE (from average wave)

figure('DefaultAxesFontSize',18); hold on
    title([num2str(VSDI.ref) '-cond' num2str( condition) ])

    for roi_i = makeRow(selroi)

    %AVERAGE
    roiwaves = percFwave(:,roi_i,:); % for all trials
    averoiwave = mean(squeeze(roiwaves), 2);

    if savewaves

        %plot result
        plot(timebase_adj, averoiwave, 'DisplayName',roilabels{roi_i}, 'linewidth', 2)
        ylabel('%F')
    end

    temp = devo_peak2peak(averoiwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
    percFpeak(roi_i,ci) = temp.peakminusbasel*10;
    percFwmean(roi_i,ci) = temp.wmean*10; 
    clear roiwaves averoiwave temp
end



end % METHOD 2

%%  -------------------------------------------
% METHOD 3 - Extract wave from average movie using 
% -------------------------------------------
switch measure_mehod
    case '%F_meth3_avewise_roi'

            % -------------------------------------------
    % SELECT CASES  AND AVERAGE MOVIE
    % -------------------------------------------
    sel_trials = find(VSDI.condition(:,1)==condition);

    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end

  
    %----------------------------------------------------------------
    % WAVES
    %----------------------------------------------------------------        
        
     % --------------------------------------------------------------------------
    % GET AVERAGE MOVIE. Use timerange set
    % --------------------------------------------------------------------------

    movieave = mean(movies(:,:,idxrange,sel_trials),4);
    movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
    ...(that normally corresponds to the background)

        % -------------------------------------------------------
        % CALCULATE %F WAVE FOR EACH ROI
        % -------------------------------------------------------
        
        % DEPRECATED 22/03/22 
        meanF0 = squeeze(mean(VSDmov.F0(:,:,sel_trials),3));
        
        for roii = makeRow(selroi)
            roimask = masks(:,:,roii);
            %                 allroi_waves(:,roii,ci) = roi_TSave(movieave,roimask);
            allroi_waves(:,roii) = roi_TSave_percF_roiwise(movieave,roimask, meanF0); %CHANGE TO NOT %F
        end %for roi_i
        
        
        
        % PLOT
        figure
        
        back = VSDI.backgr(:,:,VSDI.nonanidx(1));
        
        sp1 = subplot(1,2,1);
        switch roikind
            case 'circle'
                centers = VSDI.roi.circle.center(selroi, :) ;
                roicirc_preview_multiple_cmap(back, centers, VSDI.roi.circle.R, sp1, cmap);
                
            case 'anat'
                roi_preview_multiple(back, VSDI.roi.manual_poly(selroi,:), sp1);
        end
        
        sp1.Visible = 0;
        
        sp2= subplot(1,2,2);
        
        pos1 = get(sp1,'Position');
        pos2 = get(sp2,'Position');
        pos3= [pos2(1) pos2(2) pos1(3) pos1(4)];
        set(sp2, 'Position',pos3)
        
        nroi = length(selroi);
        roicolors= roi_colors();
        
        hold on
        i = 0; %counter for colors
        for roii = selroi
            i = i+1;
            waveroi = movmean(allroi_waves(wave.start:wave.end,roii),5);
            plot(VSDI.timebase(wave.start:wave.end), waveroi , 'linewidth', 1.3, 'Color', cmap(i,:), 'displayname' , VSDI.roi.labels_circ{roii});
            clear waveroi
            % legend(selroinames{:}, 'Location', 'northeast')
        end
        
    

end % METHOD 3

%% FINAL TOUCHS

        legend

        xlim([timebase_adj(1) timebase_adj(2)])
        ylabel('%F')
        
        sgtitle([num2str(VSDI.ref), ref_movie, 'rej', num2str(reject_on), '-cond:' num2str( condition)])
        
        savename= ['waves' num2str(VSDI.ref) ref_movie '_cond' num2str( condition) '_rej'  num2str(reject_on) '_' num2str(numel(sel_trials)) 'trials_' ];
        
        if savewaves
            saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
            
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
            print(fullfile(savein, [savename '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them %STILL TO TEST!
            
            close all
        end

% 
% if savewaves
%     savenamew = [num2str(VSDI.ref) '_' ref_movie '_rej' num2str(reject_on) '_cond' num2str( condition)];
%     %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
%     set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%     print(fullfile(savein, [analysisref '_' savenamew '_WAVES' '.jpg']),'-r300','-djpeg') % prints it as you see them (for poster: -r900)
%     close
% end%

