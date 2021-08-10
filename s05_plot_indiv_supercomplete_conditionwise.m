%% INDIVIDUAL TRIALS: SUPER COMPLETE (ecg + GS + waves + roi preview + connect matrix + tiles)
... and it does it condition_wise, saving each condition in a separate folder 

clear 
user_settings; 
% set(0,'DefaultFigureVisible','off')

for nfish = [11 12]
VSDI = TORus('load',nfish);
% spike = TORus('loadspike', nfish);

% Load movie for tiles
movie_ref = '_06filt3'; %@ SET input movie
VSDmov = TORus('loadmovie',nfish,movie_ref);

% load waves to plot: 
VSDroiTS =TORus('loadwave',nfish); 
roi2plotidx = [1 3 5 7 11];  %@ SET

% Tiles and waves SETTINGS 
    absthres =  0.05;%@ SET threshold to zero-out the movie...;
    ylim_waves = [-0.1 0.6];%@ SET y-axis limits to plot waves
    act_clim= [-0.3 0.3];  %@ SET coloraxis of shown colors

    %Calculate connectivity matrix to plot (all rois- just to check if it's
    %all yellow)
    fcM = Fconn_matrix (VSDroiTS.filt306.data, 1:14); %all rois included

%% get all different conditions
trial_kinds = unique(VSDI.condition(:,1));
trial_kinds = trial_kinds(~isnan(trial_kinds));

for which_trials =  makeRow(trial_kinds)
    
    % config folder in which trials from same condition will be saved
    foldername = [num2str(VSDI.ref), '_', num2str(which_trials)]; 
    parentfolder = fullfile(path.rootpath, 'plot', 'indiv_trials', 'reject2' ); %@ SET
    mkdir(parentfolder, foldername); 
    
    savein = fullfile(parentfolder, foldername); 
   
%% Get index from rejected trials
rejectidx = [];

    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];

    rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];

    rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];

    rejectidx = setdiff(rejectidx, VSDI.reject.forcein);


%%
    
for triali =  makeRow(find(VSDI.condition(:,1) == which_trials) ) %loops through the trials from the condition

    figure
    trialref = VSDI.list(triali).Name(end-7:end-4);
    trialref = VSDI.list(triali).Name(end-7:end-4);
    if ismember(triali, rejectidx)
    sgtitle([num2str(VSDI.ref) ': trial ' trialref '#' num2str(VSDI.list(triali).code) '(' num2str(VSDI.list(triali).mA) 'mA) - Rejected'])
    else
    sgtitle([num2str(VSDI.ref) ': trial ' trialref '#' num2str(VSDI.list(triali).code) '(' num2str(VSDI.list(triali).mA) 'mA)'])
    end

    %% 1. Plot ECG

%     ax1 = subplot (5,5,[1 2]);
%     
%     %get values
%     [ecg_crop,tbase] = get_ecg_crop(VSDI.spike.RO(triali), 30, spike.ecg, spike.ecg_timebase);
%     plot(tbase, ecg_crop)
%     ylim([-2.5 2.5]);
%     xline(VSDI.spike.RO(triali)+VSDI.info.Sonset/1000,'--r', 'LineWidth',1.5);

    %% 2. PLOT WAVES

    ax2 = subplot (5,5,[3 4 8 9]);
    timeseries = VSDroiTS.filt306.data;
    titletxt = '(right hemisph)';
    stim_dur = VSDI.list(triali).Sdur;

    for i = 1:length(roi2plotidx)
        shownrois(i) = VSDI.roi.labels(roi2plotidx(i));
    end

    plot(VSDI.timebase, VSDroiTS.filt306.data(:,roi2plotidx,triali), 'Linewidth', 1.5);
    ylim(ylim_waves)
    xlim([VSDI.timebase(1) VSDI.timebase(end)])
    yline(absthres, '--k', 'Linewidth',0.8) %plots tile's visualization threshold
    yline(act_clim(2), '--k', 'Linewidth',0.8) %plots color limit (above zero)

    legend([shownrois 'color thresh'])
    if VSDI.list(triali).mA ~=0 %if there is stimulus, plot
        xline(0, '--r', 'Linewidth',1.3); hold on
        xline(VSDI.info.Sdur(triali), '--r', 'Linewidth',1.3);
    end

    %% 3. PLOT GS from trial vs al the rest

    ax3 = subplot (5,5,[6 7]);
    %select trials from same condition:
    same_idx = VSDI.condition(:,1) == VSDI.condition(triali,1);
    GS = VSDroiTS.filt306.GS; %select those for plotting together
    titletxt = 'GS trial vs all trials';

    plot(VSDI.timebase, GS(:,same_idx), 'Linewidth', 0.1); hold on %'color', [1 0.8 0.3]
    plot(VSDI.timebase, GS(:,triali), '-k', 'Linewidth', 2);
    ylim(ylim_waves)
    xlim([VSDI.timebase(1) VSDI.timebase(end)])
    yline(absthres, '--k', 'Linewidth',0.8)
    if VSDI.list(triali).mA ~=0 %if there is stimulus, plot
        xline(0, '--r', 'Linewidth',1.3); hold on
        xline(VSDI.info.Sdur(triali), '--r', 'Linewidth',1.3);
    end
    xlim([VSDI.timebase(1) VSDI.timebase(end)])
%     annotation('textbox',[.17 .68 .1 .1],'String','GS','FitBoxToText','on');

%% 4. PLOT roi preview

    ax4 = subplot (5,5,5);
    
    imagesc(VSDI.backgr(:,:,triali)); colormap('bone'); hold on
    roicolors= roi_colors();

    for nroi = 1:size(VSDI.roi.manual_poly,1)
        coord = VSDI.roi.manual_poly{nroi};
        plot(coord(:,1), coord(:,2), 'color', roicolors(nroi,:), 'LineWidth', 1); hold on
    end

 %% 5. PLOT CONNECTIVITY MATRIX

    ax5 = subplot (5,5,10);
    timeser = VSDroiTS.filt306.data(:,:,triali);
    imagesc(fcM(:,:,triali));
    set(gca, 'clim', [-1 1])
    colormap(ax5, polarmap(parula))
    
  %% 6. Plot Tiles 
    movie = VSDmov.data(:,:,1:end-1,triali); %movie to plot

        % ... what and how to plot
        t2lpot = linspace(-10, 600, 12);

        idx2plot = find_closest_timeidx(t2lpot, VSDI.timebase);
        realt2plot = VSDI.timebase(idx2plot);
        % ...noise threshold settings
    %     baseline = 1:100; %range of frames to consider baseline (idx)
    %     SDfactor = 4;
        %...absolute threshold settings
        method = 'absvalue';


    % ABSOLUTE THRESHOLD 
    [movie_thres, alphachan] = movie_absthresh(movie,method,absthres); 
    % NOTE_DEV: test how it has to be the alphachannel

    plotsites = [9:20]; 

    for ii= 1:length(idx2plot)
        ploti = plotsites(ii);
        % OVERLAID PLOT
        background = VSDI.backgr(:,:,triali);
        % frame = movie(:,:,frame2plot);
        frame = movie_thres(:,:,idx2plot(ii));
        framealpha = alphachan(:,:,idx2plot(ii));

        ax(ploti) = subplot(5,4,ploti);
        title([num2str(realt2plot(ii)) 'ms']);
        plot_framesoverlaid(frame,background, framealpha ,0, ax(ploti), act_clim, 1); 
        % imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar
        ax(ploti).XTick = []; ax(ploti).YTick = [];
    end

    trialref = num2str(VSDI.trialref(triali)); trialref = [trialref(end-2:end) 'A']; 
    name2save = fullfile(savein,['indiv' num2str(VSDI.ref) '_trial' trialref '.jpg']);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12 12]);
    saveas(gcf,name2save,'jpg')
% 
    close
clear GS realt2plot idx2plot t2lpot which_trials

end %trial
end % which_trials
end %nfish
set(0,'DefaultFigureVisible','on')
