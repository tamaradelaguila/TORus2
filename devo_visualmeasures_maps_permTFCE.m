%% s06_analysis1: FRAME-WISE STATISTICAL COMPARISON (TFCE CORRECTION)
clear

nfish =  12%@ SET

% clear
user_settings

[VSDI] = TORus('load',nfish);

temp = TORus('loadmovie',nfish,'_06filt3');
movies = temp.data(:,:,1:end-1,:);

%% GET CONDITION CODES
%             cond_codes = unique(VSDI.condition(:,1));
%             cond_codes=  cond_codes(~isnan(cond_codes));
%             cond_codes = setdiff(cond_codes, 0); %delete code=0 if present
cond_codes =[401:404]; %for nfish = 6 (#210412)

reject_on = [1]  %@ SET

setting.manual_reject = 1; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET


out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_measurements' ; %@ SET

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

%% SELECT CASES
sel_trials= [];

for condi = cond_codes %to make sure that only the conditions of interest are computed)
    condtrials = makeCol(find(VSDI.condition(:,1)==condi));
    sel_trials  = [sel_trials; condtrials];
end
sel_trials= sort(sel_trials);

if reject_on  %@ SET
    sel_trials = setdiff(sel_trials, rejectidx);
end

%% 1. PEAK-2-PEAK MEASUREMENTS OF EACH TRIAL MOVIE


window.min = [-100 100];
window.max = [0 600];
window.movsum = 50;
window.basel = [-100 0];

method = 'movsum';

%--------------------------------------
% APPLY FUNCTION TO EACH TRIAL AND PLOT
%--------------------------------------
tic

for triali = makeRow(sel_trials) % we loop through conditions because we need to substract the control block, and also we ensure only needed trials are computed
    
    movtrial = squeeze(movies(:,:,:,triali));
    
    for rowi = 1:size(movtrial,1)
        for coli = 1:size(movtrial,2)
            wave = squeeze(movtrial(rowi, coli, :));
            output = devo_peak2peak(wave, VSDI.timebase, window,[], method, 0);
            
            peak2peak(rowi,coli,triali) = output.p2p_value;
            peakminusbasel(rowi,coli,triali) = output.peakminusbasel;
            peaklat(rowi,coli,triali) = output.peaklat_ms;
            p2plat(rowi,coli,triali) = output.p2plat_ms;
            onset30_latency_ms(rowi,coli,triali) = output.onset30_latency_ms;
            onsetnoise_ms(rowi,coli,triali) = output.onsetnoise_ms;
            noisethresh(rowi,coli,triali) = output.noisethresh;
            
            
                    
            waveW = wave(output.peakidx(1):output.peakidx(2));
            waveslope = diff(waveW);
            meanslope = mean(waveslope);
                    
            meanslope(rowi, coli,triali) = meanslope;

            clear output meanslope wave 
            
        end %coli
    end %rowi
    
    display(triali)
    
end %triali

t2 = toc
%         blob()


%2. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
%         maxval.peak2peak= max(abs(frames.peak2peak(:)));
%         maxval.peakminusbasel= max(abs(frames.peakminusbasel(:)));
%         maxval.onsetnoise_ms= max(abs(frames.onsetnoise_ms(:)));
%         maxval.noisethresh= max(abs(frames.noisethresh(:)));
%
%         c_lim.peak2peak = [-maxval.peak2peak maxval.peak2peak];
%         c_lim.peakminusbasel = [-maxval.peakminusbasel maxval.peakminusbasel];
%         c_lim.onsetnoise_ms = [0 maxval.onsetnoise_ms];
%         c_lim.noisethresh = [0 maxval.noisethresh];
%
%         BVmap = colormap_loadBV();


%% PERMUTATION TFCE of the measures from all trials

%--------------------------------------
% PERMUTATION IN 'peakminusbasel'
%--------------------------------------

j = 1;
% load(fullfile(cd,'data','grouplist.mat'));

for condi = makeRow(cond_codes)
    
    condA = condi;
    condB = force0ending(condi);
    
    %% SELECT CONDITIONS TO COMPARE AND LOAD MOVIES
    
    [idxA] = find(VSDI.condition(:,1)==condA);
    [idxB] = find(VSDI.condition(:,1)==condB);
    
    if reject_on
        idxA = setdiff(idxA, rejectidx);
        idxB = setdiff(idxB, rejectidx);
    end
    
    % cond_def = ['code',num2str(condA),'-', num2str(VSDI.list(idxA(1)).mA),'mA''minus-control' ];
    cond_def = ['code',num2str(condA),'-', num2str(condB) ];
    
    %--------------------------------------
    %
    %--------------------------------------
    
    %peakminusbasel
%     framesA = peakminusbasel(:,:,idxA);
%     framesB =  peakminusbasel(:,:,idxB);
    
    
        framesA = peakminusbasel(:,:,idxA);
    framesB =  peakminusbasel(:,:,idxB);

    %% SETTINGS FOR PLOTTING AND SAVING
    % for saving plots
    pathplot = fullfile(path.rootpath ,'plot','test_measurements','TFCEperm', 'peakminusbase_perm'); %@ SET
    savename = fullfile(pathplot,['peakminusbasel_perm', grouplist{nfish}(5:end-4),cond_def,'.jpg']); %@ SET
    
    load(fullfile(path.data, 'grouplist'))
    
    % For the permutation test:
    % nchoosek(24, 12)
    n_perm = 10000;  %@ SET
    
    % For Plotting:
    % act_clim= [-4 4]; %@SET coloraxis of the shown colors
    p_thresh = 0.001;  %@SET theshold p-value for showing diffmap values
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% COMPUTATION:  PERMUTATION + TFCE - Diffmap with p-value based threshold
    
    % nA = length(idxA); nB = length(idxB);
    
    Data{1} = permute( framesA(:,:,:), [3 1 2]);
    Data{2} = permute( framesB(:,:,:), [3 1 2]);% control/baseline condition. Trials have to be in the first dimension
    
    % APPLY FUNCTION: PERMUTATION + TFCE
    results = ept_TFCE(Data, 'i', n_perm); % independent trials
    
    % Maps to plot
    maps.peakminusbase.diffmap(:,:,j) = squeeze(mean(Data{1}) - mean(Data{2})); %the first dimension is trials
    maps.peakminusbase.Tobs(:,:,j) = results.Obs;
    maps.peakminusbase.Pmap(:,:,j) =results.P_Values;
    maps.peakminusbase.alphamap(:,:,j) = maps.Pmap(:,:,j) < p_thresh;
    
    clear Data results frame
        
    if reject_on
        name= [num2str(VSDI.ref) 'PEAK-B: cond' num2str(condA) 'minus' num2str(condB) 'n' num2str(n_perm) 'TFCEpermut' 'p' num2str(p_thresh),'(reject)'];
    else
        name= [num2str(VSDI.ref) 'PEAK-B: cond' num2str(condA) 'minus' num2str(condB) 'n' num2str(n_perm) 'TFCEpermut' 'p' num2str(p_thresh)];
    end
    
    j = j+1;
end

% pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/frame_perm';
% save(fullfile(pathsave,name),'maps')

%% PLOT STATISTICAL DIFFERENCE BETWEEN CONDITIONS (P-THRESHOLDED)

% for ploti = 1:length(cond_codes)
%     subplot (3,3,ploti)
%     imagesc(maps.diffmap(:,:,ploti))
%     idxcond = find(VSDI.condition(:,1) == cond_codes(ploti)); idxcond = idxcond(1);
%     title(['c=',num2str(cond_codes(ploti)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'])
%     %     colormap()
%     axis image
% end

for ploti = 1:length(cond_codes)
    clims= [0 0.05];
    maps = maps.peakminusbase;
    ax(ploti) = subplot(3, 3, ploti);
    logicalpha = maps.Pmap(:,:,ploti)< 0.01;
    %    imagesc(imtiles(:,:,ploti))
    test_plot_framesoverlaid(maps.Pmap(:,:,ploti),VSDI.backgr(:,:,VSDI.nonanidx(1)), logicalpha ,0, ax(ploti), clims, 1, parula);
%         test_plot_framesoverlaid(maps.Pmap(:,:,ploti),VSDI.backgr(:,:,VSDI.nonanidx(1)), maps.alphamap(:,:,ploti) ,0, ax(ploti), clims, 1, parula);

    
    idxcond = find(VSDI.condition(:,1) == cond_codes(ploti)); idxcond = idxcond(1);
    tit= ['c=',num2str(cond_codes(ploti)), '(', num2str(VSDI.condition(idxcond,4)), 'mA)'];
    
    set(gca,'XColor', 'none','YColor','none')
    ax(ploti).Title.String = tit; 
%     colormap(parula)
    
end


ax(9) = subplot(3,3,9)
imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)))
colormap(ax(9), bone)
axis image
% 
% title_info = ['diffmap' , name]; %@ CHANGE!!!
% %SAVE
% saveas(gcf, fullfile(pathsave,[title_info '.jpg']), 'jpg')
% close all
% clear maps  name pathsave moviesA moviesB idxA idxB nA nB cond_def


% end %fish loop

%% Created: 19/02/21
% FROM: Gent2 code 'pipel09_framewise1', updated: 17/11/20
