%% R-ONSET as described by (Reynaud 2011)

% IT NEED EXTRA SMOOTHING TO WORK

%% ON ROI-TIMESERIES
clear
user_settings
nfish = 3;%@ SET
VSDI = TORus('load',nfish);

refmovie = '_07filt4'; %@ SET
VSDroiTS = TORus('loadwave',nfish, refmovie);

timeser = VSDroiTS.filt407.data;

idx_trial = VSDI.nonanidx;
% idx_trial = find(VSDI.condition(:,1)==303);

idx_Stim = find_closest_timeidx(0, VSDI.timebase); 

diffseries = diff(timeser, [], 1);
nthresh = mean(diffseries(1:idx_Stim,:,:)) + 2.5*std(diffseries(1:idx_Stim,:,:)); 
nthresh = squeeze(nthresh);

for triali = makeRow(VSDI.nonanidx)
    for nroi = 1:length(VSDI.roi.labels)
        temp = diffseries(idx_Stim:end,nroi,triali); %get only post stimulus values
        temp_thresh = nthresh(nroi,triali);
        idx_onset = find(temp>temp_thresh , 1, 'first');
        idx_onset = idx_onset + idx_Stim; %in temp it was cropped from the moment of Stimulus arrival, so we have to sum it up to get the index
        if isempty(idx_onset); idx_onset = length(VSDI.timebase-1);end
        Ronset(nroi,triali) = VSDI.timebase(idx_onset);
    end
end

for triali = makeRow(VSDI.nonanidx)
    for nroi = 1:length(VSDI.roi.labels)

        plot(VSDI.timebase,timeser(:,nroi,triali)); xline(Ronset(nroi, triali))
title([num2str(triali) , VSDI.roi.labels{nroi}])
pause
    end
end
    
    
    
%% ON AVERAGED MOVIES (on single movies seems to be uneffective)

user_settings
nfish = 3;%@ SET
VSDI = TORus('load',nfish);


refmovie = '_07filt4'; %@ SET
VSDmov = TORus('loadmovie',nfish, refmovie);


ntrial = length(VSDI.list);

idx_trial = find(VSDI.condition(:,1)==303);

idx_onset = find_closest_timeidx(0, VSDI.timebase); 

% DERIVATIVE FOR EACH POINT
for triali = idx_trial
    
    movie= VSDmov.data(:,:,1:end-1,idx_trial);

    diffmovie(:,:,:,triali)= diff(movie, [], 3);

%     thresh = mean(diffmovie(:,:,1:idx_onset),3) + 2*std(diffpre, [],3);

end

% threshold matrix for each point
thresh = mean(diffmovie(:,:,:,:),3) + 2.5* std(diffmovie(:,:,:,:),[],3);

%preview
for ii= 50:150
   imagesc(diffmovie_ave(:,:,ii)); 
   set(gca, 'clim',[0 0.05])
   title(num2str(ii))
   pause
end

% first time above threshold
for triali = idx_trial'
    for rowi = 1:size(VSDI.backgr,1)
        for coli = 1:size(VSDI.backgr,2)
            timeser = squeeze(VSDmov.data(rowi,coli,1:end-1,triali));
            onsetidx = find( timeser> thresh(rowi,coli) ,1, 'first');
            if isempty(onsetidx); onsetidx= length(VSDI.timebase); end
            Onset(rowi, coli, triali) = VSDI.timebase(onsetidx);
        end 
    end
end

imagesc(mean(Onset,3)); colorbar; set(cga, 'clim',[0 350]); 
% max value




% THRESHOLD mean + 2.5SD

% R onset: time to reach the threshold


