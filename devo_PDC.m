%% PARTIAL DIRECTED COHERENCE + BIOGRAPH
% see papers: Pascucci 2020, Baccalá 2000
% https://github.com/PscDavid/dynet_toolbox
clear
user_settings

free
nfish = 12; 
[VSDI] =TORus('load',nfish);
[VSDroiTS]= TORus('loadwave',nfish);
tsdata = VSDroiTS.circ_filt306.data;

% select roi to analyze
    % Rroi = [1 3 5 7 9 11 13]; %[1 3 5 7 11];  %@ SET; :
    % Lroi= [2 4 6 8 10 12 14];
    bothroi = [1 2 3 4 5 6 7 8 11 12 13 14]; %excluding dld-c
    
rois = bothroi;
    
reject_on= 1;

cond_codes =[400:404];
cond_exp = 404;
cond_control= 400;

%% SETTINGS
% REJECTION INDEX

rejectidx = [];
    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
    rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];
    rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];
    rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
    
rejectidx = sort(unique(rejectidx));

% SELECT ONLY CONDITION CASES
sel_trials= [];

for condi = cond_codes %to make sure that only the conditions of interest are computed)
    condtrials = makeCol(find(VSDI.condition(:,1)==condi));
    sel_trials  = [sel_trials; condtrials];
end

sel_trials= sort(sel_trials);

if reject_on  %@ SET
    sel_trials = setdiff(sel_trials, rejectidx);
end

    
%% PDI COMPARISON SETTINGS 
idxA =intersect(find(VSDI.condition(:,1) ==cond_exp) , sel_trials); %take the condition from the selected trials (already cleaned)
idxB =intersect(find(VSDI.condition(:,1) ==cond_control) , sel_trials); %take the condition from the selected trials (already cleaned)
cond_def= strcat(num2str(VSDI.condition(idxA(1),1)), 'minus', num2str(VSDI.condition(idxB(1),1)));

% get labels
j = 1;
for i = rois
labels{j} = VSDI.roi.labels{i};
j = j+1;
end

% STOK-filter settings
popt = 6;
srate = 1000/VSDI.timeabs(2);
% frange= 1:round(srate/2);frange = frange';
frange= 1:40; frange= frange';

len = length(VSDI.timebase);

%condition exp
Yc1 = permute(tsdata(1:len,rois,idxA), [3,2,1]);

SKc1 = dynet_SSM_STOK(Yc1,popt);
sk_PDCc1 = dynet_ar2pdc(SKc1,srate,frange,'sPDC',[],[],1);
% dynet_connplot(sk_PDCc1,VSDI.timebase,frange,[],[],[],[],0)
dynet_connplot(sk_PDCc1,VSDI.timebase,frange,labels,[],[],[],0)

%condition 0
Yc0 = permute(tsdata(1:len,rois,idxB), [3,2,1]);

SKc0 = dynet_SSM_STOK(Yc0,popt);
sk_PDCc0 = dynet_ar2pdc(SKc0,srate,frange,'sPDC',[],[],1);
dynet_connplot(sk_PDCc0,VSDI.timebase,frange,[],[],[],[],0)
% sgtitle(['STOK']) % only with >MATLAB R2019a

% %% PDC from all fish
% for nfish = 2:7
%     
% [VSDI] =TORus('load',nfish);
% [VSDroiTS]= TORus('loadwave',nfish);
% tsdata = VSDroiTS.circ_filt306.data;
% % Input to the algorithm
%  
% idxA =intersect(find(VSDI.condition(:,1) ==cond_exp) , sel_trials); %take the condition from the selected trials (already cleaned)
% idxB =intersect(find(VSDI.condition(:,1) ==cond_control) , sel_trials); %take the condition from the selected trials (already cleaned)
% cond_def= strcat(num2str(VSDI.condition(idxA(1),1)), 'minus', num2str(VSDI.condition(idxB(1),1)));
% 
% % STOK-filter settings
% popt = 6;
% srate = 1000/VSDI.timeabs(2);
% frange= 1:round(srate/2);frange = frange';
% % frange= 1:10; frange= frange';
% 
% len = length(VSDI.timebase);
% 
% %condition exp
% Yc1 = permute(tsdata(1:len,rois,idxA), [3,2,1]);
% 
% SKc1 = dynet_SSM_STOK(Yc1,popt);
% sk_PDCc1 = dynet_ar2pdc(SKc1,srate,frange,'sPDC',[],[],1);
% dynet_connplot(sk_PDCc1,VSDI.timebase,frange,[],[],[],[],0)
% 
% %condition 0
% Yc0 = permute(tsdata(1:len,rois,idxB), [3,2,1]);
% 
% SKc0 = dynet_SSM_STOK(Yc0,popt);
% sk_PDCc0 = dynet_ar2pdc(SKc0,srate,frange,'sPDC',[],[],1);
% dynet_connplot(sk_PDCc0,VSDI.timebase,frange,[],[],[],[],0)
% 
% % store in all-fish structures
% PDCall_c1(nfish,:,:,:,:)= sk_PDCc1; %fish*roi*roi*freq*time
% PDCall_c0(nfish,:,:,:,:)= sk_PDCc0; %fish*roi*roi*freq*time
% 
% clearvars -except PDCall_c1 PDCall_c1
% end

%% PERMUTATION TEST - TFCE correction for the maps
% clear
% 
% save_results = 0; 
% plotnsave_results = 1;
% pathplot = fullfile(rootpath, 'plots', 'stats', 'PDC'); %SET
% 
% for n = 1:length()
%     for m = 1:10
%         if n~=m
% 
%  % DATA TO COMPARE (between groups or conditions)
% %SET
% % Data has to be "Trials x Channels(pixels) x Samples(time)"
% Data{1} = squeeze(PDCall_c1(:,n,m,:,:)) ; %fish*roi*roi*freq*time
% Data{2} = squeeze(PDCall_c0(:,n,m,:,:)); % control/baseline condition. Trials have to be in the first dimension
% 
% % APPLY FUNCTION
% results = ept_TFCE(Data, 'i', 5000); close
% close
% 
%     if save_results
%     TFCEresults_PDC{n,m} = results.P_Values ; %SET
%     save(fullfile(rootpath,'data_structures', 'stats', ['TSrasterTFCE', struct_list{fish}(5:end)]), 'TFCEresults_p10');
%     end
% 
% % Maps to plot
% Diffmap = squeeze(mean(Data{1}) - mean(Data{2})); %the first dimension is trials
% Tobs = results.Obs;
% Pmap =results.P_Values;
% % Plot Settings
% difflim = [-.15 .15];
% tlim = [-4 4];
% load(fullfile(rootpath,'functions','colormap_BV.mat'))
% roicmap= roi_colors();
% 
%     if plotnsave_results
%     fig1= figure;
%     supertitle= strcat(num2str(VSDI.ref),VSDI.roi.labels{roi_idx}, '-', cond_def);
%     sgtitle(supertitle);
% 
%     ax1 = subplot(2,2,1);
%     imagesc(VSDI.timebase, [],Diffmap);
%     colormap(ax1, polarmap(BVmap))
%     colorbar
%     set(gca,'clim',difflim)
%     xlabel('time(ms)'); ylabel('clust: med to lat')
%     title('Diffmap -no thresh (%\Delta F/F_0)')
% 
%     ax2 = subplot(2,2,2);
%     imagesc(VSDI.timebase, [],Tobs);
%     colormap(ax2, polarmap())
%     colorbar
%     set(gca,'clim',tlim)
%     xlabel('time(ms)'); ylabel('clust: med to lat')
%     title('T-obs (non-corrected)')
% 
%     ax3 = subplot(2,2,3);
%     plot_logPmap(Pmap,0, VSDI.times(1),VSDI.design.Sonset, [], [],ax3);
%     close
%     xlabel('time(ms)'); ylabel('clust: med to lat')
%     title('TFCE P-values (purple)')
%     
%     ax4 = subplot(2,2,4);
%     imagesc(VSDI.backgr(:,:,1)); hold on
%     colormap(ax4, 'bone')
%     f = fill(VSDI.roi.manual_poly{roi_idx,1}(:,1), VSDI.roi.manual_poly{roi_idx,1}(:,2), roicmap(roi_idx,:));
%     alpha(f, 0.5); hold off
%     
%     %SET
%     name = fullfile(pathplot,['p10_TFCE', struct_list{fish}(5:end-4),'roi',num2str(roi_idx),'-',cond_def,'.jpg']);
%     saveas(gcf, name, 'jpg')
%     close all
%     end
%  end %nroi
%         end %if statement (to avoid computing the diagonal)
%         end %m
% end %n

%% GET ADJACENCY MATRIX FROM PDC for a certain time window and frequency window
% 
% freq_range = 1:83; %freq range to plot
% time_ms = [0 50]; %@ SET
% 
% t1= dsearchn(time_ms(1)); t2= dsearchn(time_ms(end));
% 
% timeidx_range= t1:t2;
% 
% PDCmat = sk_PDCc1; 
% 
% [ADJmat] = PDCtoADJ(PDCmat, timeidx_range, freq_range);

%% DYNAMIC CONNECTIVITY: MOVING WINDOW

% To have comparable graph edges, the moving window has to be applied in
% two steps: one to get the adjacency matrices and 

% Step 1: get 3D adjacency matrix
PDCmat = sk_PDCc1; 

freq_range = 1:83; 
window = 5; halfw = floor(window/2);
time_range_ms = [-100 200];

    t1= dsearchn(VSDI.timebase, time_range_ms(1)); 
    t2= dsearchn(VSDI.timebase, time_range_ms(end)); 

time_range = [t1,t2];

j=1; %will be each step from the sliding window
for i = time_range(1)+halfw : time_range(2)-halfw
    timeref(j)= VSDI.timebase(i);
    timeidx_range = i-halfw:i+halfw;
    ADJmat(:,:,j) = PDCtoADJ(PDCmat, timeidx_range, freq_range);
    j=1+j;
end

% Step 2: create graph for each and save into movie format
nframes = size(ADJmat,3);
thr = 0.25; 
thrADJmat = ADJmat; thrADJmat(thrADJmat<thr) = 0; % Threshold matrix
[max_w, min_w] = get_maxnonzeromin(thrADJmat); %max and non-zero minimum

backgr = squeeze(VSDI.backgr(:,:,VSDI.nonanidx(1)));

% get labels
j = 1;
for i = rois
labels{j} = VSDI.roi.labels{i};
j = j+1;
end

coor = VSDroiTS.circ_filt306.circleroi.center(rois,:);

for fr = 1:nframes
temp = thrADJmat(:,:,fr)/10; %scaled so the weight are similar to a correlation matrix; otherwise the edges will be too thick
gtitle = strcat('t=',num2str(timeref(fr)));

fig = plot_graphTOR(temp,coor,max_w,min_w, thr, backgr, labels, gtitle);

biograph_frames(fr)=getframe(fig);
close
end

% Watch frames step by step
for fr = 1:nframes
imagesc(biograph_frames(fr).cdata)
pause
end 

% Play movie
h = axes;
% h.Position = [0 -0.5 1 1];
movie(h,f,1,20, [-.5 -.5 .5 .5])

%% Last Update: 09/06/21