%% s06_analysis1: FRAME-WISE STATISTICAL COMPARISON (TFCE CORRECTION) 

%% LOAD DATA
clear

user_settings
% for nfish= 6
nfish = 6;
refmovie = '_06filt3';

load(fullfile(path.data, 'grouplist.mat'))
VSDI = TORus('load',nfish);
VSDmov = TORus('loadmovie',nfish, refmovie);

exp_cond = [2001 2002 2003]; 
control = [2000];

%% SELECT CONDITIONS TO COMPARE AND LOAD MOVIES
for rowi = 1 %nº of rows in the condition to loop
for condA =2001  exp_cond(rowi,:) %@ SET
% condB = control(rowi); %@ SET
condB = control; %@ SET

[idxA] = find(VSDI.condition(:,1)==condA); 
[idxB] = find(VSDI.condition(:,1)==condB); 

% cond_def = ['code',num2str(condA),'-', num2str(VSDI.list(idxA(1)).mA),'mA''minus-control' ];
cond_def = ['code',num2str(condA),'-', num2str(condB) ];

moviesA = VSDmov.data(:,:,:,idxA);
moviesB = VSDmov.data(:,:,:,idxB);

%% SETTINGS FOR PLOTTING AND SAVING
% for saving plots
pathplot = fullfile(path.rootpath ,'plot', 'frame_perm'); %@ SET
savename = fullfile(pathplot,['movie3', grouplist{nfish}(5:end-4),cond_def,'.jpg']); %@ SET

load(fullfile(path.data, 'grouplist'))
savename = fullfile(pathplot,['movie3_diffmap', grouplist{nfish}(end-5:end),cond_def,'.jpg']); %@ SET


% Time-points to analyze
chosen_tpoints = linspace(-10, 800, 16) ; %@SET 3 timepoints (ms) from timebase to analyse

% For the permutation test:
% nchoosek(24, 12)
n_perm = 10000;  %@ SET

% For Plotting:
% act_clim= [-4 4]; %@SET coloraxis of the shown colors
p_thresh = 0.05;  %@SET theshold p-value for showing diffmap values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert time 'ms' into index
[~,t_idx] = min(abs(VSDI.timebase - chosen_tpoints));
timesample = VSDI.timebase(t_idx); % actual times-samples that will be used

%% COMPUTATION: 16 PERMUTATION + TFCE - Diffmap with p-value based threshold

nA = length(idxA); nB = length(idxB);

for ifr=  1:16
frame = t_idx(ifr);
Data{1} = permute( squeeze(moviesA(:,:,frame,:)), [3 1 2]); 
Data{2} = permute( squeeze(moviesB(:,:,frame,:)), [3 1 2]);% control/baseline condition. Trials have to be in the first dimension

% APPLY FUNCTION: PERMUTATION + TFCE
results = ept_TFCE(Data, 'i', n_perm); % independent trials

% Maps to plot
maps.diffmap(:,:,ifr) = squeeze(mean(Data{1}) - mean(Data{2})); %the first dimension is trials
maps.Tobs(:,:,ifr) = results.Obs;
maps.Pmap(:,:,ifr) =results.P_Values;
maps.alphamap(:,:,ifr) = maps.Pmap(:,:,ifr) < p_thresh; 

clear Data results frame
disp(ifr)
end
maps.background = VSDI.backgr(:,:,idxA(1));
maps.clim = [-0.7 0.7];
maps.timepoints = timesample;

name= [num2str(VSDI.ref) 'cond' num2str(condA) 'minus' num2str(condB) 'n' num2str(n_perm) 'TFCEpermut' 'p' num2str(p_thresh)];
pathsave = path.data;
save(fullfile(pathsave,name),'maps')

%% PLOT STATISTICAL DIFFERENCE BETWEEN CONDITIONS (P-THRESHOLDED) FOR THE 16 FRAMES
title_info = ['diffmap' num2str(VSDI.ref), ':', cond_def,'.p=', num2str(p_thresh)]; %@ CHANGE!!!
plot_16perm_results('diffmap_alpha', maps, title_info)
    %SAVE
        saveas(gcf, fullfile(pathsave,[title_info '.jpg']), 'jpg')
    close all
clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def
end % condA loop
end % rowi loop

% end %fish loop
%% Created: 19/02/21
% FROM: Gent2 code 'pipel09_framewise1', updated: 17/11/20
