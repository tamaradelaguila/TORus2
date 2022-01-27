%% s02 FILTERING (further preprocessing) and crop
% Crop masks should have been drawn in the s01_importmov_basic_preproces

clear
user_settings;
nfish = 2;

VSDI = TORus('load',nfish);

%% FILT2:
%
% %@ SET FILTERS PARAMETERS :
% tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
% mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
% meanpix = 3;
% medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
%
%
% % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% % complete list)
% inputRef =  '_02diff';
% outputRef = '_05filt2'; %@ SET
%
% % Load input movie
% [inputStruct] = TORus('loadmovie',nfish,inputRef);
% inputdata=inputStruct.data;
%
%
% % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
%
% % Preallocate in NaN
% filtmov1 = NaN(size(inputdata));
% filtmov2 = NaN(size(inputdata));
% filtmov3 = NaN(size(inputdata));
% filtmov4 = NaN(size(inputdata));
%
%
% % 2.1. Tcnst
% for triali = makeRow(VSDI.nonanidx)
% tempmov = inputdata(:,:,:,triali);
% filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
% clear tempmov
% end
%
% % 2.2. Spatial Filter (mean)
% for triali = makeRow(VSDI.nonanidx)
% tempmov = filtmov1(:,:,:,triali);
% filtmov2(:,:,:,triali) = filter_spatial(tempmov, meanpix);
% clear tempmov
% end
%
% % 2.3. Median spatial filter
% for triali = makeRow(VSDI.nonanidx)
% tempmov = filtmov2(:,:,:,triali);
% filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
% clear tempmov
% end
%
% % 2.4. Cubic filter
% for triali = makeRow(VSDI.nonanidx)
% tempmov = filtmov3(:,:,:,triali);
% filtmov4(:,:,:,triali) = filter_cubicBV(tempmov);
% clear tempmov
% end
%
% % % 2.5. Crop background
% % for triali = makeRow(VSDI.nonanidx)
% %     tempmov = filt4(:,:,:,triali);
% %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
% %     clear tempmov
% % end
%
% % SET definitive filtered movie that will be stored (do not forget to set the filters
% % accordingly
% filtmov_def = filtmov4;
%
% % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% % structure used to apply new changes in
% VSDmov.ref = inputStruct.ref;
% VSDmov.movieref= outputRef;
% VSDmov.data = filtmov4;
% VSDmov.times = inputStruct.times;
% %@ SET !!! according to the filters applie (append as many as needeed)
% VSDmov.hist = inputStruct.hist;
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'tcnst = 10'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'spatialmean = 3'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'median = 3'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'cubicBV'; %@ SET
% % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
% TORus('savemovie', VSDmov, VSDmov.movieref);
%
%
% %% FILT1:
% %@ SET FILTERS PARAMETERS :
% tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
% gauss_kernel = 1; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
% mpix = 3; % pixels size for MEDIAN SPATIAL FILTER
%
% % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% % complete list)
% inputRef =  '_02diff';
% outputRef = '_04filt1'; %@ SET
%
% % Load input movie
% [inputStruct] = TORus('loadmovie',nfish,inputRef);
% inputdata=inputStruct.data;
%
% % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
%
% % Preallocate in NaN
% filtmov1 = NaN(size(inputdata));
% filtmov2 = NaN(size(inputdata));
% filtmov3 = NaN(size(inputdata));
% filtmov4 = NaN(size(inputdata));
%
% % 2.1. Tcnst
% for triali = makeRow(VSDI.nonanidx)
% tempmov = inputdata(:,:,:,triali);
% filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,10);
% clear tempmov
% end
%
% % 2.2. Gaussian spatial (not into function yet)
% for triali = makeRow(VSDI.nonanidx)
%     tempmov = filtmov1(:,:,:,triali);
%      filtmov2(:,:,:,triali) = filter_gauss2D(tempmov, gauss_kernel);
%      clear tempmov
% end
%
% % 2.3. Median spatial filter
% for triali = makeRow(VSDI.nonanidx)
%     tempmov = filtmov2(:,:,:,triali);
%      filtmov3(:,:,:,triali) = filter_median(tempmov, mpix);
%      clear tempmov
% end
%
% % 2.4. Crop background
% for triali = makeRow(VSDI.nonanidx)
%     tempmov = squeeze(filtmov3(:,:,:,triali));
%     filtmov4(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
%     clear tempmov
% end
%
% % SET definitive filtered movie that will be stored (do not forget to set the filters
% % accordingly
% filtmov_def = filtmov3;
%
% % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% % structure used to apply new changes in
% VSDmov.ref = inputStruct.ref;
% VSDmov.movieref= outputRef;
% VSDmov.data = filtmov_def;
% VSDmov.times = inputStruct.times;
% %@ SET !!! according to the filters applie (append as many as needeed)
% VSDmov.hist = inputStruct.hist;
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'tcnst = 10'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'gauss = 1'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'median = 3'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
% TORus('savemovie', VSDmov, VSDmov.movieref);

%% FILT3:

%@ SET FILTERS PARAMETERS :
tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
% mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
meanpix = 9;
medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER


% 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_03diff_perc';
outputRef = '_06filt3'; %@ SET

% Load input movie
[inputStruct] = TORus('loadmovie',nfish,inputRef);
inputdata=inputStruct.data;


% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

% Preallocate in NaN
filtmov1 = NaN(size(inputdata));
filtmov2 = NaN(size(inputdata));
filtmov3 = NaN(size(inputdata));
filtmov4 = NaN(size(inputdata));


% 2.1. Tcnst
for triali = makeRow(VSDI.nonanidx)
    tempmov = inputdata(:,:,:,triali);
    filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
    clear tempmov
    disp(triali)
end

% 2.2. Spatial Filter (mean)
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov1(:,:,:,triali);
    filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
    clear tempmov
end

% 2.3. Median spatial filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov2(:,:,:,triali);
    filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
    clear tempmov
end

% 2.4. Cubic filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov3(:,:,:,triali);
    filtmov4(:,:,:,triali) = filter_cubicBV(tempmov);
    clear tempmov
end



% % 2.5. Crop background
% for triali = makeRow(VSDI.nonanidx)
%     tempmov = filt4(:,:,:,triali);
%     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
%     clear tempmov
% end

% SET definitive filtered movie that will be stored (do not forget to set the filters
% accordingly
filtmov_def = filtmov4;

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.movieref= outputRef;
VSDmov.data = filtmov4;
VSDmov.times = inputStruct.times;
%@ SET !!! according to the filters applie (append as many as needeed)
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = 'tcnst = 10'; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = 'spatialmean = 9'; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = 'median = 3'; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = 'cubicBV'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
TORus('savemovie', VSDmov, VSDmov.movieref);


blob()
% %% FILT4:
% clear
% user_settings;
%
% for nfish = 4
%
% VSDI = TORus('load',nfish);
% %@ SET FILTERS PARAMETERS :
% tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
% % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
% meanpix = 3;
% medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
%
%
% % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% % complete list)
% inputRef =  '_03diff_perc';
% outputRef = '_07filt4'; %@ SET
%
% % Load input movie
% [inputStruct] = TORus('loadmovie',nfish,inputRef);
% inputdata=inputStruct.data;
%
%
% % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
%
% % Preallocate in NaN
% filtmov1 = NaN(size(inputdata));
% filtmov2 = NaN(size(inputdata));
% filtmov3 = NaN(size(inputdata));
% filtmov4 = NaN(size(inputdata));
%
%
% % 2.1. Tcnst
% for triali = makeRow(VSDI.nonanidx)
% tempmov = inputdata(:,:,:,triali);
% filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
% clear tempmov
% disp(triali)
% end
%
% % 2.2. Spatial Filter (mean)
% for triali = makeRow(VSDI.nonanidx)
% tempmov = filtmov1(:,:,:,triali);
% filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
% clear tempmov
% end
%
% % 2.3. Median spatial filter
% for triali = makeRow(VSDI.nonanidx)
% tempmov = filtmov2(:,:,:,triali);
% filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
% clear tempmov
% end
%
% % 2.4. Cubic filter
% for triali = makeRow(VSDI.nonanidx)
% tempmov = filtmov3(:,:,:,triali);
% filtmov4(:,:,:,triali) = filter_cubicBV(tempmov);
% clear tempmov
% end
%
%
%
% % % 2.5. Crop background
% % for triali = makeRow(VSDI.nonanidx)
% %     tempmov = filt4(:,:,:,triali);
% %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
% %     clear tempmov
% % end
%
% % SET definitive filtered movie that will be stored (do not forget to set the filters
% % accordingly
% filtmov_def = filtmov4;
%
% % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% % structure used to apply new changes in
% VSDmov.ref = inputStruct.ref;
% VSDmov.movieref= outputRef;
% VSDmov.data = filtmov4;
% VSDmov.times = inputStruct.times;
% %@ SET !!! according to the filters applie (append as many as needeed)
% VSDmov.hist = inputStruct.hist;
% VSDmov.hist{length(VSDmov.hist)+1,1} = ['tcnst = ' num2str(tcnst) ]; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean =' num2str(meanpix)]; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = ['median = ' num2str(medianpix)]; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'cubicBV'; %@ SET
% % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
% TORus('savemovie', VSDmov, VSDmov.movieref);
%
% clear VSDmov filtmov1 filtmov2 filtmov3 filtmov4 filtmov_def
% end

%% 09FILT3 :
clear
user_settings;

for nfish = [5:6]
    user_settings;
    
    VSDI = TORus('load',nfish);
    
    %@ SET FILTERS PARAMETERS :
    tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
    % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
    meanpix = 9;
    medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
    
    
    % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
    % complete list)
    inputRef =  '_08diff_perc_f0pre';
    outputRef = '_09filt3'; %@ SET
    
    % Load input movie
    [inputStruct] = TORus('loadmovie',nfish,inputRef);
    inputdata=inputStruct.data;
    
    
    % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
    
    % Preallocate in NaN
    filtmov1 = NaN(size(inputdata));
    filtmov2 = NaN(size(inputdata));
    filtmov3 = NaN(size(inputdata));
    filtmov4 = NaN(size(inputdata));
    
    
    % 2.1. Tcnst
    for triali = makeRow(VSDI.nonanidx)
        tempmov = inputdata(:,:,:,triali);
        filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
        clear tempmov
        disp(triali)
    end
    
    % 2.2. Spatial Filter (mean)
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov1(:,:,:,triali);
        filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
        clear tempmov
    end
    
    % 2.3. Median spatial filter
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov2(:,:,:,triali);
        filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
        clear tempmov
    end
    
    % 2.4. Cubic filter
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov3(:,:,:,triali);
        filtmov4(:,:,:,triali) = filter_cubicBV(tempmov);
        clear tempmov
    end
    
    
    
    % % 2.5. Crop background
    % for triali = makeRow(VSDI.nonanidx)
    %     tempmov = filt4(:,:,:,triali);
    %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    %     clear tempmov
    % end
    
    % SET definitive filtered movie that will be stored (do not forget to set the filters
    % accordingly
    filtmov_def = filtmov4;
    
    % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
    % structure used to apply new changes in
    VSDmov.ref = inputStruct.ref;
    VSDmov.movieref= outputRef;
    VSDmov.data = filtmov4;
    VSDmov.times = inputStruct.times;
    %@ SET !!! according to the filters applie (append as many as needeed)
    VSDmov.hist = inputStruct.hist;
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'tcnst = 10'; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'spatialmean = 9'; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'median = 3'; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'cubicBV'; %@ SET
    % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
    TORus('savemovie', VSDmov, VSDmov.movieref);
    
    blob()
    clear
end

%% 11 FILT2 :
clear
user_settings;

for nfish =  [11:12]
    user_settings;
    
    VSDI = TORus('load',nfish);
    
    %@ SET FILTERS PARAMETERS :
    tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
    % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
    meanpix = 3;
    medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
    
    
    % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
    % complete list)
    inputRef =  '_10diff_perc_f0pre_crop';
    outputRef = '_11filt4'; %@ SET
    
    % Load input movie
    [inputStruct] = TORus('loadmovie',nfish,inputRef);
    inputdata=inputStruct.data;
    
    
    % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
    
    % Preallocate in NaN
    filtmov1 = NaN(size(inputdata));
    filtmov2 = NaN(size(inputdata));
    filtmov3 = NaN(size(inputdata));
    filtmov4 = NaN(size(inputdata));
    
    
    % 2.1. Tcnst
    for triali = makeRow(VSDI.nonanidx)
        tempmov = inputdata(:,:,:,triali);
        filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
        clear tempmov
        disp(triali)
    end
    
    % 2.2. Spatial Filter (mean)
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov1(:,:,:,triali);
        filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
        clear tempmov
    end
    
    % 2.3. Median spatial filter
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov2(:,:,:,triali);
        filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
        clear tempmov
    end
    
    
    % % 2.5. Crop background
    % for triali = makeRow(VSDI.nonanidx)
    %     tempmov = filt4(:,:,:,triali);
    %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    %     clear tempmov
    % end
    
    % SET definitive filtered movie that will be stored (do not forget to set the filters
    % accordingly
    filtmov_def = filtmov3;
    
    % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
    % structure used to apply new changes in
    VSDmov.ref = inputStruct.ref;
    VSDmov.movieref= outputRef;
    VSDmov.data = filtmov_def;
    VSDmov.times = inputStruct.times;
    %@ SET !!! according to the filters applie (append as many as needeed)
    VSDmov.hist = inputStruct.hist;
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'tcnst = 10'; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'spatialmean = 3'; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'median = 3'; %@ SET
    % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
    TORus('savemovie', VSDmov, VSDmov.movieref);
    
    blob()
    clear
end

%% 12 FILT5 :
clear
user_settings;

for nfish =  [4 9]
    user_settings;
    
    VSDI = TORus('load',nfish);
    
    %@ SET FILTERS PARAMETERS :
    tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
    % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
    meanpix = 9;
    medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
    
    
    % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
    % complete list)
    inputRef =  '_10diff_perc_f0pre_crop';
    outputRef = '_12filt5'; %@ SET
    
    % Load input movie
    [inputStruct] = TORus('loadmovie',nfish,inputRef);
    inputdata=inputStruct.data;
    
    
    % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
    
    % Preallocate in NaN
    filtmov1 = NaN(size(inputdata));
    filtmov2 = NaN(size(inputdata));
    filtmov3 = NaN(size(inputdata));
    filtmov4 = NaN(size(inputdata));
    
    
    % 2.1. Tcnst
    for triali = makeRow(VSDI.nonanidx)
        tempmov = inputdata(:,:,:,triali);
        filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
        clear tempmov
        disp(triali)
    end
    
    % 2.2. Spatial Filter (mean)
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov1(:,:,:,triali);
        filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
        clear tempmov
    end
    
    % 2.3. Median spatial filter
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov2(:,:,:,triali);
        filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
        clear tempmov
    end
    
    
    % % 2.5. Crop background
    % for triali = makeRow(VSDI.nonanidx)
    %     tempmov = filt4(:,:,:,triali);
    %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    %     clear tempmov
    % end
    
    % SET definitive filtered movie that will be stored (do not forget to set the filters
    % accordingly
    filtmov_def = filtmov3;
    
    % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
    % structure used to apply new changes in
    VSDmov.ref = inputStruct.ref;
    VSDmov.movieref= outputRef;
    VSDmov.data = filtmov_def;
    VSDmov.times = inputStruct.times;
    %@ SET !!! according to the filters applie (append as many as needeed)
    VSDmov.hist = inputStruct.hist;
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'tcnst = 10'; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'spatialmean = 3'; %@ SET *** IT WAS WRONG: it was 9!!! ***
    VSDmov.hist{length(VSDmov.hist)+1,1} = 'median = 3'; %@ SET
    % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
    TORus('savemovie', VSDmov, VSDmov.movieref);
    
    blob()
    clear
end


%% 13 (FILT 6) - BASIC BRAINVISION ON DIFFperc MOVIES, POSTERIOR CROP
% ----------------------------------------

clear
user_settings;
nfish = 2;

VSDI = TORus('load',nfish);

%@ SET FILTERS PARAMETERS :
meanpix = 9;% Time-constant for TEMPORAL SMOOTHING
medianpix = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING

% 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_03diff_perc'; %it was diffperc
outputRef = '_13filt6'; %@ SET

% Load input movie
[inputStruct] = TORus('loadmovie',nfish,inputRef);
inputdata=inputStruct.data;

% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

% Preallocate in NaN
filtmov1 = NaN(size(inputdata));
filtmov2 = NaN(size(inputdata));
filtmov3 = NaN(size(inputdata));
filtmov4 = NaN(size(inputdata));


% 2.2. Spatial Filter (mean)
for triali = makeRow(VSDI.nonanidx)
    tempmov = inputdata(:,:,:,triali);
    filtmov1(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
    clear tempmov
end

% 2.3. Median spatial filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov1(:,:,:,triali);
    filtmov2(:,:,:,triali) = filter_median(tempmov, medianpix);
    clear tempmov
end

% 2.4. Cubic filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov2(:,:,:,triali);
    filtmov3(:,:,:,triali) = filter_cubicBV(tempmov);
    clear tempmov
end

% 2.5. Crop background
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov3(:,:,:,triali);
    filtmov4(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    clear tempmov
end

% SET definitive filtered movie that will be stored (do not forget to set the filters
% accordingly
filtmov_def = filtmov4;

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.movieref= outputRef;
VSDmov.data = filtmov_def;
VSDmov.times = inputStruct.times;
%@ SET !!! according to the filters applie (append as many as needeed)
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean = ' num2str(meanpix)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = ['median = ' num2str(medianpix)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = 'cubicBV'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
TORus('savemovie', VSDmov, VSDmov.movieref);
blob()

%% 14 (FILT 6) - BASIC BRAINVISION ON DIFF MOVIES, POSTERIOR CROP
% ----------------------------------------

clear
user_settings;
nfish = 11;

VSDI = TORus('load',nfish);

%@ SET FILTERS PARAMETERS :
meanpix = 9;% Time-constant for TEMPORAL SMOOTHING
medianpix = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING

% 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_02diff';
outputRef = '_14filt6'; %@ SET

% Load input movie
[inputStruct] = TORus('loadmovie',nfish,inputRef);
inputdata=inputStruct.data;

% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

% Preallocate in NaN
filtmov1 = NaN(size(inputdata));
filtmov2 = NaN(size(inputdata));
filtmov3 = NaN(size(inputdata));
filtmov4 = NaN(size(inputdata));


% 2.2. Spatial Filter (mean)
for triali = makeRow(VSDI.nonanidx)
    tempmov = inputdata(:,:,:,triali);
    filtmov1(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
    clear tempmov
end

% 2.3. Median spatial filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov1(:,:,:,triali);
    filtmov2(:,:,:,triali) = filter_median(tempmov, medianpix);
    clear tempmov
end

% 2.4. Cubic filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov2(:,:,:,triali);
    filtmov3(:,:,:,triali) = filter_cubicBV(tempmov);
    clear tempmov
end

% 2.5. Crop background
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov3(:,:,:,triali);
    filtmov4(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    clear tempmov
end

% SET definitive filtered movie that will be stored (do not forget to set the filters
% accordingly
filtmov_def = filtmov4;

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.movieref= outputRef;
VSDmov.data = filtmov_def;
VSDmov.times = inputStruct.times;
%@ SET !!! according to the filters applie (append as many as needeed)
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean = ' num2str(meanpix)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = ['median = ' num2str(medianpix)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = 'cubicBV'; %@ SET
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
TORus('savemovie', VSDmov, VSDmov.movieref);
blob()

%% 12 FILT5 :
clear
user_settings;

for nfish =  [2 8 11 12]
    user_settings;
    
    VSDI = TORus('load',nfish);
    
    %@ SET FILTERS PARAMETERS :
    tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
    % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
    meanpix = 9;
    medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
    
    
    % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
    % complete list)
    inputRef =  '_02diff';
    outputRef = '_15filt5'; %@ SET
    
    % Load input movie
    [inputStruct] = TORus('loadmovie',nfish,inputRef);
    inputdata=inputStruct.data;
    
    
    % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
    
    % Preallocate in NaN
    filtmov1 = NaN(size(inputdata));
    filtmov2 = NaN(size(inputdata));
    filtmov3 = NaN(size(inputdata));
    filtmov4 = NaN(size(inputdata));
    
    
    % 2.1. Tcnst
    for triali = makeRow(VSDI.nonanidx)
        tempmov = inputdata(:,:,:,triali);
        filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
        clear tempmov
        disp(triali)
    end
    
    % 2.2. Spatial Filter (mean)
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov1(:,:,:,triali);
        filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
        clear tempmov
    end
    
    % 2.3. Median spatial filter
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov2(:,:,:,triali);
        filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
        clear tempmov
    end
    
    
    % % 2.5. Crop background
    % for triali = makeRow(VSDI.nonanidx)
    %     tempmov = filt4(:,:,:,triali);
    %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    %     clear tempmov
    % end
    
    % SET definitive filtered movie that will be stored (do not forget to set the filters
    % accordingly
    filtmov_def = filtmov3;
    
    % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
    % structure used to apply new changes in
    VSDmov.ref = inputStruct.ref;
    VSDmov.movieref= outputRef;
    VSDmov.data = filtmov_def;
    VSDmov.times = inputStruct.times;
    %@ SET !!! according to the filters applie (append as many as needeed)
    VSDmov.hist = inputStruct.hist;
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['tcnst =' num2str(tcnst)]; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean = ' num2str(meanpix)]; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['median =' num2str(medianpix)]; %@ SET
    % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
    TORus('savemovie', VSDmov, VSDmov.movieref);
    
    blob()
    clear
end

%% _17filt5 (from diff_f0preS-10frames)
clear
user_settings;

for nfish =  [8:12]
    user_settings;
    
    VSDI = TORus('load',nfish);
    
    %@ SET FILTERS PARAMETERS :
    tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
    % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
    meanpix = 9;
    medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
    
    
    % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
    % complete list)
    inputRef =  '_16diff_f0pre';
    outputRef = '_18filt6'; %@ SET
    
    % Load input movie
    [inputStruct] = TORus('loadmovie',nfish,inputRef);
    inputdata=inputStruct.data;
    
    
    % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
    
    % Preallocate in NaN
    filtmov1 = NaN(size(inputdata));
    filtmov2 = NaN(size(inputdata));
    filtmov3 = NaN(size(inputdata));
    filtmov4 = NaN(size(inputdata));
    
    
    % 2.1. Tcnst
    for triali = makeRow(VSDI.nonanidx)
        tempmov = inputdata(:,:,:,triali);
        filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
        clear tempmov
        disp(triali)
    end
    
    % 2.2. Spatial Filter (mean)
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov1(:,:,:,triali);
        filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
        clear tempmov
    end
    
    % 2.3. Median spatial filter
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov2(:,:,:,triali);
        filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
        clear tempmov
    end
    
    
    % % 2.5. Crop background
    % for triali = makeRow(VSDI.nonanidx)
    %     tempmov = filt4(:,:,:,triali);
    %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    %     clear tempmov
    % end
    
    % SET definitive filtered movie that will be stored (do not forget to set the filters
    % accordingly
    filtmov_def = filtmov3;
    
    % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
    % structure used to apply new changes in
    VSDmov.ref = inputStruct.ref;
    VSDmov.movieref= outputRef;
    VSDmov.data = filtmov_def;
    VSDmov.times = inputStruct.times;
    %@ SET !!! according to the filters applie (append as many as needeed)
    VSDmov.hist = inputStruct.hist;
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['tcnst =' num2str(tcnst)]; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean = ' num2str(meanpix)]; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['median =' num2str(medianpix)]; %@ SET
    % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
    VSDmov.F0 = inputStruct.F0; %@ SET
    
    TORus('savemovie', VSDmov, VSDmov.movieref);
    
    blob()
    clear
end

%% _18filt6 (from diff_f0preS-10frames) + bleach removal
clear
user_settings;

for nfish =  [2]
    user_settings;
    
    VSDI = TORus('load',nfish);
    
    %@ SET FILTERS PARAMETERS :
    tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
    % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
    meanpix = 9;
    medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
    cutoff = 0.1 ; %Hz for the butterworth high-pass BLEACH-LIKE FILTER 
    
    % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
    % complete list)
    inputRef =  '_16diff_f0pre';
    outputRef = '_18filt6'; %@ SET
    
    % Load input movie
    [inputStruct] = TORus('loadmovie',nfish,inputRef);
    inputdata=inputStruct.data;
    
    
    % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
    
    % Preallocate in NaN
    filtmov1 = NaN(size(inputdata));
    filtmov2 = NaN(size(inputdata));
    filtmov3 = NaN(size(inputdata));
    filtmov4 = NaN(size(inputdata));
    
    
    % 2.1. Tcnst
    for triali = makeRow(VSDI.nonanidx)
        tempmov = inputdata(:,:,:,triali);
        filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
        clear tempmov
        disp(triali)
    end
    
    % 2.2. Spatial Filter (mean)
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov1(:,:,:,triali);
        filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
        clear tempmov
    end
    
    % 2.3. Median spatial filter
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov2(:,:,:,triali);
        filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
        clear tempmov
    end
    
    % 2.4. Slow drift removal
    for triali = makeRow(VSDI.nonanidx)
        tempmov = filtmov3(:,:,:,triali);
        filtmov4(:,:,:,triali) = filter_bleach_butterhigh(tempmov, cutoff, VSDI.info.stime);
        clear tempmov
    end
    % % 2.5. Crop background
    % for triali = makeRow(VSDI.nonanidx)
    %     tempmov = filt4(:,:,:,triali);
    %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
    %     clear tempmov
    % end
    
    % SET definitive filtered movie that will be stored (do not forget to set the filters
    % accordingly
    filtmov_def = filtmov4;
    
    % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
    % structure used to apply new changes in
    VSDmov.ref = inputStruct.ref;
    VSDmov.movieref= outputRef;
    VSDmov.data = filtmov_def;
    VSDmov.times = inputStruct.times;
    %@ SET !!! according to the filters applie (append as many as needeed)
    VSDmov.hist = inputStruct.hist;
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['tcnst =' num2str(tcnst)]; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean = ' num2str(meanpix)]; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['median =' num2str(medianpix)]; %@ SET
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['butterhigh=' num2str(cutoff) 'Hz']; %@ SET
    VSDmov.F0 = inputStruct.F0; %@ SET
    
    % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
    TORus('savemovie', VSDmov, VSDmov.movieref);
    
    blob()
    clear
end
blob(); pause(0.1); blob()

%% _19filt7 (from diff_f0preS-10frames) only bleach removal
clear
user_settings;

for nfish =  [2 8 9 10 11 12]
    user_settings;
    
    VSDI = TORus('load',nfish);
    
    %@ SET FILTERS PARAMETERS :
    cutoff = 0.1 ; %Hz for the butterworth high-pass BLEACH-LIKE FILTER 
    
    % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
    % complete list)
    inputRef =  '_16diff_f0pre';
    outputRef = '_19filt7'; %@ SET
    
    % Load input movie
    [inputStruct] = TORus('loadmovie',nfish,inputRef);
    inputdata=inputStruct.data;
    
    
    % 2. PERFORM COMPUTATIONS:
    
    % Preallocate in NaN
    filtmov1 = NaN(size(inputdata));
    
    % 2.4. Slow drift removal
    for triali = makeRow(VSDI.nonanidx)
        tempmov = inputdata(:,:,:,triali);
        filtmov1(:,:,:,triali) = filter_bleach_butterhigh(tempmov, cutoff, VSDI.info.stime);
        clear tempmov
    end
    
    % SET definitive filtered movie that will be stored (do not forget to set the filters
    % accordingly
    filtmov_def = filtmov1;
    
    % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
    % structure used to apply new changes in
    VSDmov.ref = inputStruct.ref;
    VSDmov.movieref= outputRef;
    VSDmov.data = filtmov_def;
    VSDmov.times = inputStruct.times;
    %@ SET !!! according to the filters applie (append as many as needeed)
    VSDmov.hist = inputStruct.hist;
    VSDmov.hist{length(VSDmov.hist)+1,1} = ['butterhigh=' num2str(cutoff) 'Hz']; %@ SET
    VSDmov.F0 = inputStruct.F0; %@ SET
    
    % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
    TORus('savemovie', VSDmov, VSDmov.movieref);
    
    blob()
    clear
end
blob(); pause(0.1); blob()