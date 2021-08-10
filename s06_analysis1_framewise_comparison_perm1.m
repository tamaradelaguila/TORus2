%% s06_analysis1: FRAME-WISE STATISTICAL COMPARISON (TFCE CORRECTION) 
%% OLD CODE: there is another script that allows to compute + plot 16frames
%% LOAD DATA
 clear

user_settings
nfish= 3;

load(fullfile(path.data, 'grouplist.mat'))
refmovie = '_06filt3';
VSDI = TORus('load',nfish);
VSDmov = TORus('loadmovie',nfish, refmovie);

%% SELECT CONDITIONS TO COMPARE AND LOAD MOVIES
condA = 301; %@ SET
condB = 300; %@ SET

[c_idxA] = choose_condidx(VSDI.condition(:,1), condA); 
[c_idxB] = choose_condidx(VSDI.condition(:,1), condB); 

cond_def = [num2str(condA) 'minus-control'];

moviesA = VSD.mov(:,:,:,c_idxA);
moviesB = VSD.mov(:,:,:,c_idxB);

%% SETTINGS FOR PLOTTING AND SAVING
% for saving plots
pathplot = fullfile(path.rootpath ,'plot', 'frame_perm'); %@ SET
savename = fullfile(pathplot,['movie6', grouplist{nfish}(5:end-4),cond_def,'.jpg']); %@ SET

% SET NOISE-THRESHOLDING PARAMETERS: 
% SDcrit = 1.5; %@SET noise will be considered the signal within SDcrit * ST of the baseline
% baseline_ms = [-100,0];  %@ SET: timeframe considered to be baseline 
% 
% Time-points to analyze
chosen_tpoints = [-10, 300, 600] ; %@SET 3 timepoints (ms) from timebase to analyse

% For the permutation test: 
n_perm = 5000;  %@ SET

% For Plotting:
act_clim= [-0.3 0.3]; %@SET coloraxis of the shown colors
p_thresh = 0.05;  %@SET theshold p-value for showing diffmap values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert time 'ms' into index
[~,t_idx] = min(abs(VSDI.timebase - chosen_tpoints));
timesample = VSDI.timebase(t_idx); % actual times-samples that will be used

[~,base1] = min(abs(VSDI.timebase - baseline_ms(1)));
[~,base2] = min(abs(VSDI.timebase - baseline_ms(end)));

baselidx = base1:base2;

%% COMPUTATION  
% NOISE-THRESHOLDED TRIAL-AVERAGED MOVIES
[movaveA_th, imthresA] = movie_noisethresh(mean(movieA,4),baselidx, SDcrit, 1 ) ; 
[movaveB_th, imthresB] =  movie_noisethresh(mean(movieB,4),baselidx, SDcrit, 1 ) ;

% FRAMES FROM THE THRESHOLDED MOVIE
aveframeA(:,:,1:3) = movaveA_th(:,:,t_idx);
aveframeB(:,:,1:3) = movaveB_th(:,:,t_idx);
alphaA = aveframeA~=0 ;% define alphadata for the color (activity) map: 1=show; 0=hide; 
alphaB = aveframeB~=0;

%% COMPUTATION: PERMUTATION + TFCE - Diffmap with p-value based threshold

nA = length(idxA); nB = length(idxB);

for ifr= 1:3
frame = t_idx(ifr);
Data{1} = permute( squeeze(movieA(:,:,frame,:)), [3 1 2]); 
Data{2} = permute( squeeze(movieB(:,:,frame,:)), [3 1 2]);% control/baseline condition. Trials have to be in the first dimension

% APPLY FUNCTION: PERMUTATION + TFCE
results = ept_TFCE(Data, 'i', n_perm); % independent trials

% Maps to plot
Diffmap(:,:,ifr) = squeeze(mean(Data{1}) - mean(Data{2})); %the first dimension is trials
Tobs(:,:,ifr) = results.Obs;
Pmap(:,:,ifr) =results.P_Values;
clear Data results frame
end

alphaP = Pmap<p_thresh; % Transparency for the statistical comparisons

%% PLOT STATISTICAL DIFFERENCE BETWEEN CONDITIONS (P-THRESHOLDED) FOR THE 3 FRAMES

fig = figure;
% CONDITION TONE-2000Hz
sgtitle([num2str(VSDI.ref),':', cond_defA '(up).', cond_defB, '(middle). Difference (down) [\DeltaF/F_0]'], 'Fontsize',10)

ax1 = subplot(3,3,1);
% title(['t=', num2str(timesample(1)), 'ms'], 'Fontsize', 8);
plot_framesoverlaid(movaveA_th(:,:,t_idx(1)), VSDI.backgr(:,:,1), alphaA(:,:,1) ,0, ax1, act_clim, 0); 
% imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar
ax1.XTick = []; ax1.YTick = [];

ax2 = subplot(3,3,2);
% title(['t=', num2str(timesample(2)), 'ms'], 'Fontsize', 8);
plot_framesoverlaid(movaveA_th(:,:,t_idx(2)), VSDI.backgr(:,:,1), alphaA(:,:,2), 0, ax2, act_clim, 0); 
ax2.XTick = []; ax2.YTick = [];

ax3 = subplot(3,3,3);
% title(['t=', num2str(timesample(3)),'ms'], 'Fontsize', 8) 
plot_framesoverlaid(movaveA_th(:,:,t_idx(3)), VSDI.backgr(:,:,1), alphaA(:,:,3), 0, ax3, act_clim, 1); 
ax3.XTick = []; ax3.YTick = [];

% CONDITION NO-TONE
ax4 = subplot(3,3,4);
title(['t=', num2str(timesample(1)), 'ms'], 'Fontsize', 8);
plot_framesoverlaid(movaveB_th(:,:,t_idx(1)), VSDI.backgr(:,:,1), alphaB(:,:,1) ,0, ax4, act_clim, 0); 
% imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar
ax4.XTick = []; ax4.YTick = [];

ax5 = subplot(3,3,5);
title(['t=', num2str(timesample(2)), 'ms'], 'Fontsize', 8);
plot_framesoverlaid(movaveB_th(:,:,t_idx(2)), VSDI.backgr(:,:,1), alphaB(:,:,2), 0, ax5, act_clim, 0); 
ax5.XTick = []; ax5.YTick = [];

ax6 = subplot(3,3,6);
title(['t=', num2str(timesample(3)),'ms'], 'Fontsize', 8) 
plot_framesoverlaid(movaveB_th(:,:,t_idx(3)), VSDI.backgr(:,:,1), alphaB(:,:,3), 0, ax6, act_clim, 1); 
ax6.XTick = []; ax6.YTick = [];


% PLOT STATISTICAL RESULTS (from comparing the activity in a timepoint in
% between both conditions) 

% CONDITION NO-TONE
ax7 = subplot(3,3,7);
title(['DIFFMAPS: t=', num2str(timesample(1)), 'ms'], 'Fontsize', 8);
plot_framesoverlaid(Diffmap(:,:,1), VSDI.backgr(:,:,1), alphaP(:,:,1) ,0, ax7, act_clim, 0); 
% imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar
ax7.XTick = []; ax7.YTick = [];

ax8 = subplot(3,3,8);
title(['t=', num2str(timesample(2)), 'ms'], 'Fontsize', 8);
plot_framesoverlaid(Diffmap(:,:,2), VSDI.backgr(:,:,1), alphaP(:,:,2) ,0, ax8, act_clim, 0); 
% imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar
ax8.XTick = []; ax8.YTick = [];

ax9 = subplot(3,3,9);
title(['t=', num2str(timesample(3)), 'ms. (p-threshold <' , num2str(p_thresh)], 'Fontsize', 8);
plot_framesoverlaid(Diffmap(:,:,3), VSDI.backgr(:,:,1), alphaP(:,:,3) ,0, ax9, act_clim, 1); 
% imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar
ax9.XTick = []; ax9.YTick = [];

    %SAVE
        saveas(fig, savename, 'jpg')
    close all
clear

%% Created: 19/02/21
% FROM: Gent2 code 'pipel09_framewise1', updated: 17/11/20
