% s05 EXTRACT ROI timeseries
clear

user_settings
for nfish= 11%[8 10 11 12]

VSDI = TORus('load',nfish);
VSDroiTS = TORus('loadwave',nfish);

%% CREATE STRUCTURES FROM SELECTED MOVIES

% 1. REFERENCES for input movies/output waves (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_18filt6'; 
fieldref = strcat(inputRef(4:end),inputRef(2:3)); 

%load input movie 
[inputStruct] = TORus('loadmovie',nfish,inputRef); 
movies=inputStruct.data;

% 2. PERFORM COMPUTATIONS:  EXTRACT WAVES 
nroi =length(VSDroiTS.roi.labels);
allroi_waves = NaN(length(VSDI.timebase), nroi, VSDI.nonanidx(end)); 

for triali = makeRow(VSDI.nonanidx)
    movietrial = movies(:,:,:,triali);
    for roii = 1:nroi
        roimask = VSDroiTS.roi.manual_mask(:,:,roii);
        allroi_waves(:,roii,triali) = roi_TSave(movietrial,roimask);
    end
    disp(triali)
end

% 3.SAVE in TS STRUCTURE: Save waves (VSDroiTS) structure copying some references from the movie
% structure used to apply new changes in
VSDroiTS.(fieldref).data= allroi_waves;
VSDroiTS.(fieldref).times = inputStruct.times;
VSDroiTS.(fieldref).hist = inputStruct.hist;
TORus('savewave', VSDroiTS, VSDroiTS.ref); 
clear inputStruct movies nroi allroi_waves

end

%% EXTRACT GLOBAL SIGNAL FROM CROPMASK

clear

user_settings
for nfish=  [8 10 11 12]
    
VSDI = TORus('load',nfish);
VSDroiTS = TORus('loadwave',nfish);

% 1. REFERENCES for input movies/output waves (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_06filt3'; 
fieldref = strcat(inputRef(4:end),inputRef(2:3)); 

%load input movie 
[inputStruct] = TORus('loadmovie',nfish,inputRef); 
movies=inputStruct.data;

% 2. PERFORM COMPUTATIONS:  EXTRACT WAVES 
nroi =length(VSDroiTS.roi.labels);
GS = NaN(length(VSDI.timebase), VSDI.nonanidx(end)); 
cropmask = VSDI.crop.mask;

for triali = makeRow(VSDI.nonanidx)
    movietrial = movies(:,:,:,triali);
    GS(:,triali) = roi_TSave(movietrial,cropmask);
    disp(triali)
end

% 3.SAVE in TS STRUCTURE: Save waves (VSDroiTS) structure copying some references from the movie
% structure used to apply new changes in
VSDroiTS.(fieldref).GS= cropmask;
VSDroiTS.(fieldref).GS= GS;
TORus('savewave', VSDroiTS, VSDroiTS.ref); 
clear inputStruct movies GS

end

blob()

%% s05 EXTRACT CIRCULAR ROI timeseries
clear

user_settings
for nfish=11% [8 10 11 12]

VSDI = TORus('load',nfish);
VSDroiTS = TORus('loadwave',nfish);

%% CREATE STRUCTURES FROM SELECTED MOVIES

% 1. REFERENCES for input movies/output waves (see 'z_notes.txt', point 5 for
% complete list)
inputRef = '_18filt6'; %  '_06filt3';
fieldref = strcat('circ_', inputRef(4:end),inputRef(2:3)); 

%load input movie 
[inputStruct] = TORus('loadmovie',nfish,inputRef); 
movies=inputStruct.data;

% 2. PERFORM COMPUTATIONS:  EXTRACT WAVES 
nroi =length(VSDI.roi.labels_circ);
allroi_waves = NaN(length(VSDI.timebase), nroi, VSDI.nonanidx(end)); 

for triali = makeRow(VSDI.nonanidx)
    movietrial = movies(:,:,:,triali);
    for roii = 1:nroi
        roimask = VSDI.roi.circle.mask(:,:,roii);
        allroi_waves(:,roii,triali) = roi_TSave(movietrial,roimask);
    end
    disp(triali)
end

% 3.SAVE in TS STRUCTURE: Save waves (VSDroiTS) structure copying some references from the movie
% structure used to apply new changes in
VSDroiTS.(fieldref).data= allroi_waves;
VSDroiTS.(fieldref).times = inputStruct.times;
VSDroiTS.(fieldref).hist = inputStruct.hist;
VSDroiTS.(fieldref).circleroi = VSDI.roi.circle;
VSDroiTS.(fieldref).circleroi.labels = VSDI.roi.labels_circ;

TORus('savewave', VSDroiTS, VSDroiTS.ref); 
clear inputStruct movies nroi allroi_waves

display('circular rois extracted')
blob()
end

blob(); pause(0.1); blob()

%% Updated: 30/08/21