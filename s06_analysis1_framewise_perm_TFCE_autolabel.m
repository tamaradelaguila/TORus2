%% s06_analysis1: FRAME-WISE STATISTICAL COMPARISON (TFCE CORRECTION) 

%% LOAD DATA
clear

user_settings
% for nfish= 6
nfish = 8;
refmovie = '_06filt3';

load(fullfile(path.data, 'grouplist.mat'))
VSDI = TORus('load',nfish);
VSDmov = TORus('loadmovie',nfish, refmovie);


reject_on = [ 1]  %@ SET
        
        setting.manual_reject = 0; %@ SET
        setting.GSmethod_reject = 1;  %@ SET
        setting.GSabsthres_reject = 1; %@ SET
        setting.force_include = 0; %@ SET
        
        
        out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/indiv_trials/reject2_label' ; %@ SET
        
        %% CONFIG REJECTION OPTIONS
        
        
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
        
        %% GET CONDITION CODES
            cond_codes = unique(VSDI.condition(:,1));
            cond_codes=  cond_codes(~isnan(cond_codes));
            cond_codes = setdiff(cond_codes, 0); %delete code=0 if present

            
for exp_cond =makeRow(cond_codes)
condA = exp_cond;
condB = force0ending(exp_cond);

%% SELECT CONDITIONS TO COMPARE AND LOAD MOVIES



[idxA] = find(VSDI.condition(:,1)==condA); 
[idxB] = find(VSDI.condition(:,1)==condB); 

if reject_on
    idxA = setdiff(idxA, rejectidx);
    idxB = setdiff(idxB, rejectidx);
end

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
chosen_tpoints = linspace(-10, 600, 16) ; %@SET 3 timepoints (ms) from timebase to analyse

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

if reject_on
name= [num2str(VSDI.ref) 'cond' num2str(condA) 'minus' num2str(condB) 'n' num2str(n_perm) 'TFCEpermut' 'p' num2str(p_thresh),'(reject)'];
else
name= [num2str(VSDI.ref) 'cond' num2str(condA) 'minus' num2str(condB) 'n' num2str(n_perm) 'TFCEpermut' 'p' num2str(p_thresh)];
end
pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/frame_perm';
save(fullfile(pathsave,name),'maps')

%% PLOT STATISTICAL DIFFERENCE BETWEEN CONDITIONS (P-THRESHOLDED) FOR THE 16 FRAMES
title_info = ['diffmap' , name]; %@ CHANGE!!!
plot_16perm_results('diffmap_alpha', maps, title_info)
    %SAVE
        saveas(gcf, fullfile(pathsave,[title_info '.jpg']), 'jpg')
    close all
clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def

end
% end %fish loop
%% Created: 19/02/21
% FROM: Gent2 code 'pipel09_framewise1', updated: 17/11/20
