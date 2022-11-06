% ATTENTION
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% This code is a copy from the original 'Bradycardia01_flexible_counting_code_MASTER.m' stored in a separate 'counting'
% folder. Changes in the code should be done (or updated) in the original code, and only
% parameters change should be performed in this one.

%% COMMON CODE FOR BRADYCARDIA COMPUTATION 

% The code takes list of matfiles from a csv file (source.filelist), gets the spike information from them and
% perform the computation of the bradycardic response - both by the beats counting-based method, or the %ibi -  outputting an excel
% file and a rasterplot for each kind of stimulus.

% This code can be used for multiple experimental designs:
% - Single or double stimulus designs
% - With or without RO 
% - Conditioning paradigms with or without a noisy US. 

% The code is also flexible for the IBI measure, so there are different
% ways of computing the %bradychardia: 
% - Using either a time window or a number of ibis to calculate the
% baseline(pre) or post counting 
% - Using the mean or the max value of the 'post' values, or the value of
% the ibi in a specific position
% The 'beats-based counting 

% Before using the code, the SETTINGS have to be adjusted 

% Reference for counting: 'stim' variable.
% When RO is present, the only difference is that the blank trials will be computed, but the reference
% will remain the 'stim'

% Tested for design.kind = 'single_stim'; design.RO = 0;

% Channels from spike to export into the matlab as events: spikes, stim1.
% Optional channels according to the design: stim2, RO, us

% In the csv file, the different filenames have to be writen without quotation
% marks and separated by ','. E.g.:
%   220428spike1_habit,
%   220428spike2_train,

% NOTE ABOUT BLANK TRIALS:
% both Brady and Raster for blank trials will be calculated respect to the
% estimated time where the stim arrives in stimulus-trials, so the counting
% respect to the RO arrival is always the same in all trials

% ////////////////////////////////////////////////////////////////////////////////////
% @SET:SETTINGS:
% ////////////////////////////////////////////////////////////////////////////////////

clear

source.in = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike'; 
source.out = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike/excels_brady'; 
source.filelist= 'spikelist_toprocess.csv';

design.RO = 1; %1 = with RO, 0= no RO (only spike). Will influence in whether trying to find blanks or not
design.kind = 'single_stim'; %'single_stim', 'discrim'

% RASTER SETTINGS: 'none', 'match_ibis', 'match_beats', 'custom'
% .......................................................................
config.rasterplot = 'match_beats';

% ONLY FOR config.rasterplot = 'custom':
w.raster.pre = [];
w.raster.post = [];

% ADVANCED SETTING...
% .......................................................................
w.stimdur = 0; % only used for plotting in the raster the stim offset
w.usonset = []; % ONLY FOR CONDITIONING onset respect to the CS onset for conditioning paradigms

% ADVANCED NOISE-RELATED SETTING...
% .......................................................................

% US noise (only for CONDITIONING experimental designs). 
% Set to 0 if there is no US or the US is not noisy (or if it has been cleaned manually from the spike).
preUS_noise = 0 ; % ms 
USnoise_ms = 0; % ms 

% For NOISY STIMULUS
stim_noise_ms = 0; 

% ADVANCED RO-RELATED SETTINGS (if present)
% .......................................................................
% When the RO is present, there are 2 possibilities: 
% 1 - that the RO is triggered manually from BV, in that case the RO signal in spike corresponds to the start
...of the recording and the shutdelay has to be set to 0 
% stim = RO + baseline  *(shutdelay =0)!!!
% 2 - that the RO is triggered from Spike or the stimulator, then:
% stim = RO + shutdelay + baseline

w.shutdelay = 0 ; % ms of shutter delay to count (0 if BV was ruling since the spike signals the start of the recording)
w.baseline = 0.3; % ms of baseline

% Whether we have a 'stim' event or it has to be calculated from RO (the
% latter baing config.ROreferenced = 1)
config.ROreferenced = 0; % set to 1 only if the 'stim1' vector needs to be calculated from the RO vector. Only possible for 'single_stim' design
% note: 'stim1' always needs to exist because is the reference for brady
% calculation, so if it doesnt exist by itself, it need to be computed

% /////////////////////////////////////////////////////////////////
% WINDOWS FOR COMPUTING BRADYCHARDIA
% /////////////////////////////////////////////////////////////////

% ------------------------------------------------------------------------
% WINDOWS FOR beats
% ------------------------------------------------------------------------
        w.beats.pre = 10;  % pre-stimulus or baseline time window (s)
        w.beats.post = 6;  % post-stimulus time window (s)

        % CONDITIONING DESIGNS, these windows could be used: 
%         w.beats.pre = usonset - w.preUS_noise;  % pre-stimulus or baseline time window (s) - DONT NEED TO SUBSTRACT THE preUSnoise, but just so they match
%         w.beats.post = w.usonset - w.preUS_noise;  % post-stimulus time window (s)


% ------------------------------------------------------------------------
% WINDOWS FOR COMPUTING ibis
% ------------------------------------------------------------------------

% For counting method = 'n' (number of ibi)
% ...........................................
w.ibi.pre_n= 8; % nº of spikes for the baseline. NEEDED for config.ibiproxy.pre= 'n'
w.ibi.post_n = []; % post-stimulus nº of spikes. NEEDED for config.ibiproxy.post= 'n'

% For counting method = 's' (inside window)
% ...........................................
w.ibi.pre = w.beats.pre; 
w.ibi.post =  w.beats.post; 

% For counting method = 'specific' (certain ibi)
% ...........................................
% position of that ibi (ibi=0 is the first ibi post stim)
w.ibi.locat = 1;

% ------------------------------------------------------------------------
% ADVANCED CONFIGURATION FOR %IBI CALCULATION 
% ------------------------------------------------------------------------

% (1) BASELINE IBI METHOD: 'n', 's'
%How to calculate the baseline for ibi:
% 'n' = from a determined number of spikes
% 's' - from all the spikes inside a time window (in seconds)
config.ibiproxy.pre = 'n';

% (2) POST IBI COUNTING METHOD: 'n' , 's', ''specific'
% 'n' - all 'n' number of spikes after the stimulus ATT (see below)
... senses. It needs the definition of the variable 'w.ibi.post_n'
    % - 's' - It needs the definition of the variable 'w.ibi.post'
config.ibiproxy.post = 'specific';
% ATT: 'n' will not take into account potential noise due to the
% US, so this measure only makes sense for CC with a non-noisy US or a 'single_stim' design with a
% non-noisy stim

% (3) % IBI COMPUTING METHOD: 'max', 'mean', 'specific'
% Which ibi value to take as 'post' measure:
% 'max': the maximum of the post ibi values
% 'mean' the mean of the post ibi values
% 'specific' - the % will be computed using the value of a single specific
% ibi. If ibiproxy.post = 'specific', this variable will be forced to
% 'specific'
config.ibiproxy.comp = 'specific'; % The most logical to use in  'ibiproxy.post'= 'n' or 's' is 'mean', since the baseline is a mean
% For  'ibiproxy.post'= 'specific', use 'specific' as well 

% /////////////////////////////////////////////////////////////////
% END OF CONIFGURATION 
% /////////////////////////////////////////////////////////////////

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% INTPUT CONTROL AND PROCESSING
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% ADJUST NOISE UNITS TO SECONDS
% .......................................................................
w.preUS_noise = preUS_noise/1000 ; %in seconds
w.USout= USnoise_ms/1000; % in seconds
w.Sout = stim_noise_ms/1000;  % in seconds

% ------------------------------------------------------------------------
% BUILD THE LABELS FOR LATER CHARACTERIZATION OF IBI
% ------------------------------------------------------------------------
% They will be use in the excel file, either as header or in 'params' sheet

% For the parameter definition sheet:
switch config.ibiproxy.pre
    case 'n'
        labelibi.pre = ['pre_n=', num2str(w.ibi.pre_n)];
    case 's'
        labelibi.pre = ['pre in:', num2str(w.ibi.post), 's'];
end

% For the parameter definition sheet:
switch config.ibiproxy.post
    case 'n'
        labelibi.post = ['post_n=', num2str(w.ibi.post_n) ];
    case 's'
        labelibi.post = ['post in:', num2str(w.ibi.post),'s' ];
    case 'specific'
        labelibi.post = [num2str(w.ibi.locat),'ibi post' ];
        if ~strcmpi(config.ibiproxy.comp, 'specific')
        config.ibiproxy.comp = 'specific'; %it forces the computation to be of the kind 'specific
        warning('config.ibiproxy.comp is set to "specific"')
        end
end


% For the header:
switch config.ibiproxy.comp
    case 'mean'
        labelibi.comp = 'ibi_meanpost';
    case 'max'
        labelibi.comp = 'ibi_maxpost';
    case 'specific'
        labelibi.comp = ['ibi=' num2str(w.ibi.locat)];
end

% Check that the parameters 'post' and 'comp' are compatible.
switch config.ibiproxy.comp
    case 'specific'
 if ~strcmpi(config.ibiproxy.post, 'specific')
     warning(['The ibi counting method is:' config.ibiproxy.post 'so the computing %ibi method will set to "mean"'])
     config.ibiproxy.comp = 'mean';
 end
end
% If post = 'specific', comp will be forced to be specific later in the
% code; but if post is not specific, comp cannot be ir neither. If this
% happens, 'comp' will be set to 'mean'



%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% GET INFO FROM CSV FILES
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% filelist = fullfile(source.in,strcat('spikelist',num2str(VSDI.ref),'.csv'));
% LIST OF FILES TO IMPORT AND INFO ABOUT SPIKE
filelistpath = fullfile(source.in,source.filelist);
% refpath = fullfile( source.in, [num2str(fishref) 'spikeref.csv']);

tablelist = readtable(filelistpath,'ReadVariableNames',false, 'delimiter', ',');
tablelist = table2cell(tablelist);

ref.filelist = tablelist(:,1);

% file_timeref = tablelist{:,2};


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% FROM EACH FILE: GET BRADYCHARDIA AND OUTPUT RESULTS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for fili = 1:size(ref.filelist,1)
    
    disp(['file:' num2str(fili)])
    file = ref.filelist{fili,1};
    nombre.input = file;
    nombre.output =[file(1:end-4) '_result'] ;
    
    % -------------------------------------
    % DEFINE EXCELNAME FOR EACH FILE
    % -------------------------------------
    
    excelname = fullfile(source.out, [file '_BRADYCOUNT.xls']);
    
    if isfile(excelname)
        % if file already exists:
        M = strcat ('ALREADY EXISTING FILE!!! / EL ARCHIVO YA EXISTE!!!.', newline , 'The file: "', file, '_BRADYCOUNT.xls', '" already exists in that location.' , newline, 'This code cannot override an existing file.',sprintf("\n"),  'TO SOLVE THE ERROR, change its name, its location or delete before running the code');
        error(M)
    end
    
    % -------------------------------------
    % DATA LOADING AND CONTROL OF ERRORS
    % -------------------------------------
    
    % CARGAR EL ARCHIVO Y COPIAR LOS CANALES 'spikes', 'stim' y 'us' (the
    % name os US will depend on the kind of experiment)
    
    load(fullfile(source.in, nombre.input))
    
    % VARIABLE CONTROL  (check if they exist, otherwise, the name could be
    % wrong
    
    if ~exist('spikes'); warning('check the names: event "spikes" does not seem to exist '); end
    if ~exist('US');  warning('check the names: event "US" does not seem to exist '); end
    
    if ~exist('stim1')
        if config.ROreferenced
        warning('even "stim1" does not exsit but it will be calculated from RO')
        stim1.times = RO.times +  w.shutdelay + w.baseline;
        else 
            error('event "stim1" DOES NOT EXIST: check the names of the channels before running the code; else, if you need to reference the counting respect to the RO signal, set "config.ROreferenced = 1"');
        end
    end
    
    if strcmpi(design.kind, 'discrim')
        if ~exist('stim2'); warning('check the names: event "stim2" does not seem to exist '); end
    end
    
    if design.RO
        if ~exist('RO'); error('check the names: event "RO" does not seem to exist '); end
    end
    
    
    
    % -------------------------------------
    % COPY CHANNELS AND COMPUTE
    % -------------------------------------
    
    beats= spikes.times;
    
    % GET BLANK TRIALS TIMES - those RO events without any stimulus nearby
    if design.RO
        
        % Get -provisionally- all CS to rule out from the RO-trials list that includes all
        % kinds
        switch design.kind
            case 'single_stim'
                allCS= stim1.times; % 'stim'
            case 'discrim'
                allCS = [stim1.times; stim2.times]; % 'stim'
        end
        
        stimonset = w.shutdelay + w.baseline;
        blanktrial = [];
        
        for ROi= RO.times' % store trials that have no CS
            if ~any((allCS > ROi) & (allCS< ROi + stimonset + 2))
                blanktrial = [blanktrial; ROi];
            end
        end
        
        cs.blank.event = blanktrial + stimonset;
        cs.blank.name = 'blank';
    end
    
    % Make a list of CS events to loop through (if RO, add the blank
    % trials)
    switch design.kind
        case 'single_stim'
            cs.stim1.event= stim1.times; % 'stim'
            cs.stim1.name = 'stim1';
            
            if design.RO % add to the list of stimuli the blank condition
                stimlist = {cs.stim1 cs.blank};
            else
                stimlist = {cs.stim1};
            end
            
        case 'discrim'
            cs.stim1.event= stim1.times; % 'stim'
            cs.stim2.event = stim2.times;
            
            cs.stim1.name = 'stim1';
            cs.stim2.name = 'stim2';
            
            if design.RO
                stimlist = {cs.stim1 cs.stim2 cs.blank};
            else
                stimlist = {cs.stim1 cs.stim2};
            end
    end
    
    
    try %in case that there is no US file,
        us = US.times;
    catch
        us =[];
        warning ([file  ': no US channel - therefore set to []' ])
    end
    
    clearvars -except fili source export ref w file nombre beats cs us fishref exp_code setopt stimlist rowi stim1 stim2 design config labelibi excelname
    
    % if size(stim, 1) > size (stim,2); stim = stim'; end
    
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % LOOP THROUGH CS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    rowi = 2; % first row will be for labels
    
    for CSi =  1:length(stimlist)
        currentCS = stimlist{CSi};
        stim = currentCS.event;
        
        if isempty(stim) % in case of empty 'blanks', skips the loop
            continue
        end
        % ---------------------------------------------------
        % i. CLEAN BEATS (discarding those inside noise windows)
        % ---------------------------------------------------
        
        %  1 . Get from 'beats' those inside the stimulus noise-window
        idx_discard = [];
        
        % noise from 'stim1'
        for ii = 1:length(stim1)
            temp =find(stim(ii) <beats & beats<stim(ii)+w.Sout);
            idx_discard= [idx_discard; temp] ;
            clear temp
        end
        
        
        % noise from 'stim2' (in case of CCdiscrim configuration)
        switch design.kind
            case 'discrim'
                for ii = 1:length(stim2)
                    temp =find(stim(ii) <beats & beats<stim(ii)+w.Sout);
                    idx_discard= [idx_discard; temp] ;
                    clear temp
                end
        end
        
        %  2 . Get from 'beats' those inside the US noise-window
        % If there is any US, US variable will be empty and the loop won't
        % be executed
        for jj = 1:length(us)
            temp =find(us(jj)-w.preUS_noise <beats & beats<us(jj)+w.USout+w.preUS_noise); %-0.01 is an security factor in case that
            idx_discard= [idx_discard; temp] ;
            clear temp
        end
        
        % 3. CLEAN NOISY SPIKES FROM 'beats'
        beats(idx_discard) = NaN;
        beats = beats(~isnan(beats));
        
        clear idx_discard
        
        % ---------------------------------------------------
        % GET IBI
        % ---------------------------------------------------
        ibi = diff(beats);
        %         frec_inst = 1./ibi;
        
        beats = beats(2:end); % so the index corresponds to those in the 'ibi' and 'frec_inst' vectors
        
        % ---------------------------------------------------
        % PERFORM SPIKES AND IBI COUNTING AND GET BRADYCARDIA
        % ---------------------------------------------------
        n = 0;
        for ii= 1:length(stim)
            n = n+1;
            
            % ------------------
            % 0. Flag presence of US
            % ------------------
            plus = find(stim(ii) < us & us <stim(ii)+ w.stimdur);
            if plus > 0; flag_us =1;  else flag_us = 0; end
            clear plus
            
            % ------------------
            % 1. BEATS COUNTING
            % ------------------
            
            preNwin = stim(ii) -w.beats.pre ; % ONLY PRE OPTICAL REGISTER BASELINE
            
            postNwin = stim(ii)+w.beats.post;
            
            %             idx_Npre = find(preNwin - w.ROout < beats & beats < stim(ii) - w.ROout);
            idx_Npre = find(preNwin  < beats & beats < stim(ii));
            
            idx_Npost = find(stim(ii) < beats & beats < postNwin);
            
            Npre = length(idx_Npre);
            Npost = length(idx_Npost); % in case of unequal pre and post counting windows
            Npost_adj = length(idx_Npost) * (w.beats.pre/w.beats.post); % in case of unequal pre and post counting windows
            
            percN = ((Npre - Npost) / Npre )*100;
            
            % ------------------
            % 2. IBI COUNTING
            % ------------------
            % Unlike 'beats' counting, the brady measured by '%ibi' will be
            % calculated differently according to how 'config.ibiproxy'
            % is defined
            
            % IBI PRE
            switch config.ibiproxy.pre
                case 'n'
                    idx_pre = find(beats < stim(ii),w.ibi.pre_n, 'last');
                case 's'
                    idx_pre = find(beats < stim(ii) & beats > (stim(ii) - w.ibi.pre));
                    
            end
            
            meanpre = mean(ibi(idx_pre));
            
            % IBI POST
            switch config.ibiproxy.post
                
                case 'n'
                    
                    %  OPTION 1 - all 'n' beats after the CS
                    % to Use ibi from beats after stimulus, no matter whether there
                    % is an US, w.ibi.post_n has to be defined in SETTINGS and this
                    % code ummuted. Only makes sense for 'single_stim' designs in
                    % which the stim makes no noise, or 'CC' desgins in which the
                    % US makes no noise
                    % ----------------------------
                    idx_post = find(beats > stim(jj),w.ibi.post_n, 'first');
                    meanpost = mean(ibi(idx_post));
                    maxpost = max(ibi(idx_post));
                    
                    % ----------------------------
                    
                case 's'
                    % OPTION 2 -  beats inside a window. The code has to accound for the cases in
                    % which the window is empty or have a single beat
                    % ----------------------------
                    idx_post = find(stim(ii) < beats & beats  < (stim(ii) + w.ibi.post));
                    
                    if numel(idx_post)> 0
                        meanpost = mean(ibi(idx_post));
                        maxpost = max(ibi(idx_post));
                    else
                        idx_prebeat = find( beats < stim(ii) ,1, 'last');
                        meanpost =  stim(ii) + w.stimdur - beats(idx_prebeat);
                        maxpost =  stim(ii) + w.stimdur - beats(idx_prebeat);
                        clear idx_prebeat
                    end
                    
                case 'specific' 
                    % since the first idx is ibi= 0, we take the
                    % w.ibi.locat+1 spike
                    idx_post = find(beats > stim(ii),w.ibi.locat+1, 'first'); %we first find all spikes till the one chosen
                    single_idx = idx_post(w.ibi.locat + 1); %then we keep only that one
                    singleibipost = ibi(single_idx); 
            end
            
            % ----------------------------
            
            % '%IBI'
            switch config.ibiproxy.comp
                case 'mean'
                    ibiM= meanpost;
                    perc = ((meanpost - meanpre) / meanpre)*100;
                case 'max'
                    ibiM = maxpost;
                    perc = ((maxpost - meanpre) / meanpre)*100;
                case 'specific'
                    ibiM = singleibipost;
                    perc = ((singleibipost - meanpre) / meanpre)*100;
            end
            
            % ----------------------------------------------------
            % STORE RASTER PLOT
            % ----------------------------------------------------
            
            switch config.rasterplot
                case 'match_beats'
                    % STORE FOR RASTER OF SPIKES
                    rasterRow = beats([idx_Npre' idx_Npost']) - stim(ii);
                    raster(CSi).values{n}= rasterRow;
                    raster(CSi).stim{n}= stim(ii);
                case 'match_ibis'
                    % STORE FOR RASTER OF SPIKES
                    rasterRow = beats([idx_pre' idx_post']) - stim(ii);
                    raster(CSi).values{n}= rasterRow;
                    raster(CSi).stim{n}= stim(ii);
                case 'custom'
                    % FIRST FIND THE SPIKES FROM CUSTOM WINDOW
                    preRwin = stim(ii) -w.raster.pre ; % ONLY PRE OPTICAL REGISTER BASELINE
                    postRwin = stim(ii)+w.raster.post;
                    
                    %             idx_Npre = find(preNwin - w.ROout < beats & beats < stim(ii) - w.ROout);
                    idx_Rpre = find(preRwin  < beats & beats < stim(ii));
                    idx_Rpost = find(stim(ii) < beats & beats < postRwin);
                    
                    % STORE FOR RASTER OF SPIKES
                    rasterRow = beats([idx_Rpre' idx_Rpost']) - stim(ii);
                    raster(CSi).values{n}= rasterRow;
                    raster(CSi).stim{n}= stim(ii);

            end
            
            % ------------------
            % STORE FOR EXCEL
            % ------------------
            
            %  BEATS
            % ..............................
            output{rowi,1} = round(stim(ii)); %percentage of brady from inst freq ( >0 means brady)
            output{rowi,2} = n; %number of stim in block
            output{rowi,3} = currentCS.name; %percentage of brady from inst freq ( >0 means brady)
            output{rowi,4} = flag_us; %US presence
            output{rowi,5} = 1; %for later fast inclusion/rejection
            
            output{rowi,7} = Npre; %percentage of brady from inst freq ( >0 means brady)
            output{rowi,8} = Npost; % nº of spikes in the post interval
            output{rowi,9} = round(Npost_adj,2); %nº of spikes adjusted to the 'pre' interval (rule of thumb)
            output{rowi,10} = round(percN,2); %percentage of brady from inst freq ( >0 means brady)
            
            clear preNwin postNwin idx_Npre idx_Npost  Npre Npost Npost_adj percN
            
            %  IBI
            % ..............................
            output{rowi,12} = round(meanpre,2); %percentage of brady from inst freq ( >0 means brady)
            %             output{rowi,13} = round(meanpost,2); % nº of spikes in the post interval
            output{rowi,13} = round(ibiM,2); %measure used to calculate the perc
            output{rowi,14} = round(perc,2); %percentage of brady from inst freq ( >0 means brady)
            
            rowi = rowi +1;
            clear idx_pre idx_post meanpre meanpost perc
        end %ii
        
        clear ii
        
        % ----------------------------------------------------
        % WRITE HEADER (according to the measures used in ibi)
        % ----------------------------------------------------
        
        header = {'stim(s)' 'n' 'CS' 'us' 'incl' '' 'npre' 'npost' 'npost_adj' '%beats' '' 'ibipre' labelibi.comp '%ibi'};
        
        % -------------------
        % PRINT RASTER PLOT
        % -------------------
        switch config.rasterplot
            case {'match_ibis' , 'match_beats' }
                figure
                aH = axes;
                rasterito = raster(CSi).values;
                rastervector = [];
                for n = 1:length(rasterito)
                    row = rasterito{n};
                    N = ones(size(row))*n;
                    scatter(row, N, 5, 'k' ,'marker', 'o', 'MarkerFaceColor', 'k')
                    hold on
                    
                    % values to store spikes in long format
                    l = length(row);
                    lR = length(rastervector);
                    rastervector(lR+1:lR+l, 1) = row;
                    rastervector(lR+1:lR+l, 2) = ones(size(row))*n;
                    rastervector(lR+1:lR+l, 3) =ones(size(row))*CSi;
                    rastervector(lR+1:lR+l, 4) = ones(size(row))*fili;
                    clear row
                    
                end
                
                rasterlab = num2str(round(stim));
                
                yticks(1:length(rasterito))
                yticklabels(rasterlab)
                aH.YDir = 'reverse';
                
                if strcmpi(currentCS.name , 'blank')
                xline(0, '--')
                xlabel('t(s) respect to the estimated stim arrival (if it where a non-blank trial)')
                else
                xline(0)
                xlabel('t(s) respect to the stim')
                end
                xline(w.stimdur, '--')

                
                
                ylabel('n trial')
                tit_raster = (['phase' num2str(fili) '-', currentCS.name]);
                title(tit_raster)
                
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                print(fullfile(source.out, ['RASTER-' nombre.output currentCS.name '.svg']),'-r600','-dsvg', '-painters') % prints it as you see them 
                close
                
                % ------------------------
                % OUTPUT RASTER IN  EXCEL
                % ------------------------
                % Pre-initialize with a large number of rows and as many columns as
                % trials
                L = length(raster(CSi).stim);
                l0 = 30; %random initial value
                R = NaN(L,l0);
                
                % Add the stimuli in the first
                for jj = 1:L
                    R(jj,1) = round(raster(CSi).stim{jj});
                    rasterRow = round(raster(CSi).values{jj},2);
                    lrow = length(rasterRow);
                    
                    % Before adding the row, IF the rasterRow is longer than the initial matrix (l0), make it
                    % bigger
                    if lrow+1>l0 %
                        [x,y] = size(R);
                        Rtemp = R;
                        R = NaN(L,lrow+1);
                        R(1:x, 1:y) = Rtemp;
                        l0 = lrow+1;
                    end
                    % Add the row to the matrix
                    R(jj,2:lrow+1) = rasterRow;
                    clear Rtemp rasterRow
                end
                
                writematrix(R, excelname, 'sheet', ['rast_', currentCS.name ])
                
                clear R 
                
        end
        
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     % END OF CS LOOP
    end % CSi (stim kind)
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    % -------------------------------------------
    % OUTPUT IN EXCEL
    % -------------------------------------------

    output(1,:) = header;
    
    writecell(output, excelname, 'sheet', 'bytrial')

    
    %  OUTPUT TRIALS ARRANGED BY STIM ORDER
    % ...........................................
    
    dim = size(output);
    output2(2:dim(1),:) = sortrows(output(2:end,:),1);
    output2(1,:) = output(1,:);
    writecell(output2, excelname, 'sheet', 'ordered')
    
    %  PARAMETERS
    % ...........................................
    
    labels.excel(1,1) = {date};
    labels.excel(2,1:2) = {'source'  mfilename('fullpath')};
    
    labels.excel(4,1) = {['design:', design.kind]}; 
    labels.excel(4,2) = {['RO=', num2str(design.RO)]}; 

    % STORE PARAMETERS USED IN THE OUTPUT MATRIX (PRE, POST VALUES AND STD)
    % ...........................................

    labels.excel(5:7,1) = {'parameters BTS:' ,['pre in :' num2str(w.beats.pre) 's'] , ['post in :' num2str(w.beats.post) 's']};
    
    labels.excel(9:12,1) = {'parameters IBI:' ,labelibi.pre, labelibi.post, ['method :' labelibi.comp]};
    
    labels.excel(14,1) ={['Stim_noise=', num2str(w.Sout * 1000) 'ms']}; % has to be inside the brackets, otherwise it'll yield two cells
    labels.excel(15,1) ={['USnoise=', num2str(w.USout*1000) 'ms']}; % values are multiplied to yield 'ms'
    labels.excel(16,1) ={['preUSnoise=', num2str(w.preUS_noise * 1000) 'ms']};
    
    labels.excel(18,1) ={['rasterplot:' config.rasterplot]};

    config.rasterplot = 'match_ibis'; %'match_ibis', 'match_beats', 'none'

    
    writecell(labels.excel, excelname, 'sheet', 'params')
    
    clear output output2 excelname
    
    
    
    
end %file (that is, phase, since there is one spike-file per phase)

% Update history
% 16/08/22 - include advanced configuration and IBI-calculation methods, print raster in excel sheet and ROestim option.
% 06/06/22 - updated
% 30/05/22 - Created (from other counting codes)
