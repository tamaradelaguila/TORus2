%% PARTS OF THE CODE
... [1] OUTPUTS EXCEL WITH BRAYCHARDIA AND BRAIN ACTIVITY MEASURES BOTH: INDIVIDUAL TRIALS AND MEAN GROUPED
    ... [2] PLOTS WAVES + CLOUDPLOTS + REGRESSION LINES EITHER FOR 'ALL TRIALS', FOR 'SHARKS' OR 'NO-SHARKS'
    
% TO COMPLETE THE CODE: ADD AN EXTRA LOGIC COLUMN INDICATING IN-RANGE
% CONDITION AND CALCULATE THE MEANS ACCORDING TO THAT

% The bradychardia is calculated leaving out a window with the US artifact
% (when it's tone, indicated by VSDI.condicion(:,1) == NaN, no safety window is applied)

% Ampitude measures from main rois are also reflected

% Shark/noShark classification is made such as if dldm's peak is higher than
% dm4'speak, a shark is considered. That peak is to be found in the whole triali, unlike the peaks used to measure the rois peaks
...that are to be found according to the input 'window'

clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady/groupbrady_new3_roi2')
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady/groupbrady_new3_old_roi')
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady/groupbrady_new4_old_roi') % with blank condition
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady/groupbrady_new4_roi2')
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady/groupbrady_new3_roi1_selection')

% % with blank condition

group_note = 'FINAL_redo_';  % to avoid overriding excel if a new group is loaded

selrow = [2 4 5 6]; % [2 4 5 6] SELECTED ROWS FROM THE 'groupbrady' structure

for rowi = selrow%size(groupbrady,1)
    nfish = groupbrady{rowi,1};
    
    % SET FISH AND CONDITION TO STUDY
    cond_codes = groupbrady{rowi,2};
    cond_description = groupbrady{rowi,3};

    % CHOOSE MOVIE REFERENCE ACCORDING TO THE DATAUNITS
    dataunits = '%dF'; % '%dF' 'dF'
    switch dataunits
        case 'dF'
            ref_movie = '_18filt6'; % input movie '_17filt5'
        case '%dF'
            ref_movie = '_21filt6';
    end
    
    % EXPORTING SETTINGS
    export.excel = 0; % GROUPEXCEL
    export.excel_indiv = 0;
    export.excel_singlefish = 1;
    %SETTINGS FOR PART [1] OF THE CODE ..............................
    w.pre= 10; % ventana pre-estímulo(en segundos) , de los que se va a tomar las espigas para hacer la media de frecuencia instantánea o el conteo
    w.post = w.pre/2; % sólo para el conteo
    
    outfield = 'wmean'; % 'peakminusbasel' 'wmean'
    
    
    % set(0,'DefaultFigureVisible','off')
    
    % SETTINGS FOR THE FUNCTION
    window.min = [-100 100]; %'feed-function' structure
    window.max = [0 1200]; % where to find the max peak
    window.movsum = 50; %ms
    window.basel = [-25 0]; %cambiar a  -100
    window.slope=50;
    window.wmean=[0 350];
    
    method = 'movsum';
    
    roikind = 'circle';
    
    % ATT: make sure that the 'line_roiname' is the same as the first item
    % of 'selroinames'
    selroinames =  groupbrady{rowi,5}; % the first one will be used for regression line

    % ................................................................
    
    % SETTINGS FOR PART [2] OF THE CODE: TRIALS SELECTION FOR PLOTS (do not affect the excel)...
    IBI = groupbrady{rowi,4}; % 0 , 1
    whattrials = 'no-sharks'; % 'all', 'sharks', 'no-sharks'

    % if bradyrange = [], the whole range will be used
    lim1= -5;
    step = 10;
    lim2 = 95;
    
    bradyrange = [lim1 lim2]; %bradychardia range that will be plot (and used to do the regression line)

    forceline = 3; %mimimum number of points to make the regression line: 'forceline'+1
    % ....................................................................
    
    % OUTPUT SETTINGS
    outfolder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady';
    excelname_i = fullfile(outfolder,['BRADYnew5_individ_' outfield ref_movie  '_' dataunits '_' whattrials '_' group_note '.xlsx']);
    excelname_g =  fullfile(outfolder,['BRADYnew5_grouped_' outfield ref_movie '_' dataunits '_' whattrials '_' group_note '.xlsx']);
    
    
    %% [1] COMPUTES BRADYCHARDIA COUNTING AND EXPORTS EXCEL FILE OF INDIVIDUAL TRIALS AND OF TRIALS GROUPED BY CONDITION
    
    VSDI = TORus('load',nfish);
    spike = TORus('loadspike', nfish); % ECG
    
    VSDmov = TORus('loadmovie',nfish,ref_movie);
    movies = VSDmov.data ;
    F0 = VSDmov.F0;
    %     window_sh= window; % DIFFERENT WINDOW FOR SHARK-FINDING  (the whole triali). It has to be defined inside the function because VSDI.timebase might differ among fishes
    %     window_sh.max = [0 VSDI.timebase(end-10)];
    
    %     output = {'triali', 'spiketime', 'cond','mA','', '%brady(count)','', 'IBIbasel', 'ibi0', 'ibi1','ibi2','ibi3','', 'estimation','dm4','dm4m','dm2','dldm'}; %header first
    
    %         params{8,1} = 'rejected';
    
    %         params{8,1} = rejectidx';
    excelname_1 = fullfile(outfolder, ['BRADYnew5_'  num2str(VSDI.ref)  '_'  num2str(cond_codes(end)) '_' outfield ref_movie  '_' dataunits '_' whattrials '_' group_note '.xlsx']); %excelname each

    %----------------------------------------------------------------
    % SELECT ROI
    %----------------------------------------------------------------
    switch roikind
        case 'circle'
            roi2plotidx =name2idx(selroinames, VSDI.roi.labels_circ);
            roilabels = VSDI.roi.labels_circ;
            masks =  VSDI.roi.circle.mask;
            
        case 'anat'
            roi2plotidx =name2idx(selroinames, VSDI.roi.labels);
            roilabels = VSDI.roi.labels;
            masks = VSDI.roi.manual_mask;
    end
    
    %----------------------------------------------------------------
    % SELECT triali
    %----------------------------------------------------------------
    ci = 0;
    conditrials = [];
    for condi = cond_codes
        ci = ci+1;
        
        conditrials = [conditrials; find(VSDI.condition(:,1) == condi)];
    end
    
    conditrials = sort(conditrials);
    
    %----------------------------------------------------------------
    % ADJUST TO THE TRIALS-INCLUSION CRITERIA
    %----------------------------------------------------------------

    switch whattrials
        
        case 'all'
            reject_on= 4;
            % REJECT ONLY FOR BRADYCHARDIA :
            setting.manual_reject = 0; %
            setting.GSmethod_reject = 0;  %
            setting.GSabsthres_reject = 0; %
            setting.forcein = 0; %
            setting.bradyvisual = 1;
            
            reject_brady = compute_rejectidx(VSDI, reject_on, setting);
            
            seltrials = setdiff(conditrials, reject_brady);
            
        case 'sharks'
            reject_on= 4;
            
            % REJECT ONLY FOR BRADYCHARDIA :
            setting.manual_reject = 0; %
            setting.GSmethod_reject = 0;  %
            setting.GSabsthres_reject = 0; %
            setting.forcein = 0; %
            setting.bradyvisual = 1;
            
            shark_idx  = [makeRow(VSDI.reject4.GSdeviat2sd) makeRow(VSDI.reject4.GSabs025) makeRow(VSDI.reject4.manual)];
            shark_idx = sort(unique(shark_idx));
            shark_idx=intersect(shark_idx, conditrials);
            
            reject_brady = compute_rejectidx(VSDI, reject_on, setting);
            
            seltrials = setdiff(shark_idx, reject_brady);
            
        case 'no-sharks'
            
            reject_on= 4;
            
            % REJECT ALL SHARKS :
            setting.manual_reject = 1; %
            setting.GSmethod_reject = 1;  %
            setting.GSabsthres_reject = 1; %
            setting.forcein = 0; %
            setting.bradyvisual = 1;
            
            reject_idx  = compute_rejectidx(VSDI, reject_on, setting);
            
            seltrials = setdiff(conditrials, reject_idx);
    end
    
    %
    
    ii= 2;
    for triali =  makeRow(seltrials)
        % ----------------------------------------------
        % SPIKE(ecg)-BASED MEASURES
        % ----------------------------------------------
        if  VSDI.spike.RO(triali) == 0
            flag_noROspike = 1;
        else
            flag_noROspike = 0;
        end
        
        
        if ~flag_noROspike
            
            stim = VSDI.spike.RO(triali)+VSDI.info.Sonset/1000 ; %stimulus arrival time (in spike) for this triali (RO time here has already left out the shutter delay, so we dont have to sum it up)
            
            
            % 1. LEAVE OUT STIMULUS NOISE (if it's a somatic stimulus):
            code_str = num2str( VSDI.condition(triali,1));
            % Get safety window to remove US noise
            if (VSDI.list(triali).Sdur >= 60) & (VSDI.condition(triali,4) >= 1) % EXPERIMENT-SPECIFIC PARAMETER: if it is a long train (60ms) of 1mA or more
                noisewin = 430/1000; % (430ms) empirically measured from the noise captured in the ecg
            elseif strcmpi(code_str(end), '0')
                noisewin = 0;
            else
                noisewin = VSDI.info.Sdur(triali)/1000 + 0.03; %+0.03security window (from spike, the events normally have a 0.02 minimum step from spike to spike)
            end
            
            % %         % !!!  Get safety window FOR FISH 210320 - COMMENT FOR THE REST OF FISH
            %         if (VSDI.condition(triali,1) == 303) % EXPERIMENT-SPECIFIC PARAMETER: if it is a long train (60ms) of 1mA or more
            %             noisewin = 0.41; % (430ms) empirically measured from the noise captured in the ecg
            %
            %         elseif (VSDI.condition(triali,1) == 300)
            %             noisewin = 0; % (430ms) empirically measured from the noise captured in the ecg
            %         else
            %             noisewin = 0.21; %+0.03security window (from spike, the events normally have a 0.02 minimum step from spike to spike)
            %         end
            
            
            % 2. CLEAN THE SPIKES THAT ARE DUE TO THAT TRIAL'S NOISE
            clean_spikes = spike.spikes(spike.spikes<(stim - 0.01)); %keep all before the stim (excluding another prestim safety window)
            
            if ~isnan(VSDI.condition(triali,1))% IF IT IS NaN, IT IS A TONE, so the artifact due to the stim does not exist
                clean_post = spike.spikes(spike.spikes> (stim + noisewin) );
            elseif isnan(VSDI.condition(triali,1))% If it is a tone, no 'noisewin' clearing is needed
                clean_post = spike.spikes(spike.spikes> stim);
            end
            
            clean_spikes = [clean_spikes; clean_post];
            
            % 3. GET IDX FOR BRADYCHARDIA MEASURE FROM THE CLEANED HEART
            % SPIKES
            idx_pre = find(clean_spikes > stim-w.pre  & clean_spikes < stim);
            
            idx_post =  find(clean_spikes > stim  & clean_spikes < stim+ w.post);
            idx_spike0 = find(clean_spikes> stim, 1, 'first');
            
            
            % bradycardia counting
            precount = length(idx_pre);
            postcount = length(idx_post)*2;
            brady_count = ((precount - postcount)/precount)*100;
            
            % IBI calculation
            preibi = mean(diff(clean_spikes(idx_pre)));
            
            post0 = clean_spikes(idx_spike0) - clean_spikes(idx_spike0-1);
            post1 = clean_spikes(idx_spike0+1) - clean_spikes(idx_spike0);
            post2 = clean_spikes(idx_spike0+2) - clean_spikes(idx_spike0+1);
            post3 = clean_spikes(idx_spike0+3) - clean_spikes(idx_spike0+2);
            
            % PERCENT IBI calculation
            perc0 = ((post0 - preibi) / preibi)*100;
            perc1 = ((post1- preibi) / preibi)*100;
            perc2 = ((post2 - preibi) / preibi)*100;
            perc3 = ((post3 - preibi) / preibi)*100;
            
        end % flag
        
        % ----------------------------------------------
        % BRAINVISION/WAVE-BASED MEASURES 
        % MEASURES
        % ----------------------------------------------
        movie = movies(:,:,:,triali);
        
        for nroi = makeRow(roi2plotidx)
            roimask= masks(:,:,nroi);
            
            roiwave =  roi_TSave(movie,roimask);
            temp = devo_peak2peak(roiwave, VSDI.timebase, window,[], method, 0);
            %             local_output_sh = devo_peak2peak(roiwave, VSDI.timebase, window_sh,[], method, 0);
            
            measure(nroi) = temp.(outfield);
        end
        
        % ----------------------------------------------
        % RESULS - triali-referenced so they can be subset to group
        % ----------------------------------------------
        output(triali,1) = triali;
        output(triali,2) = VSDI.trialref(triali);
        output(triali,3) = VSDI.condition(triali,1);
        output(triali,4) = VSDI.spike.RO(triali);
        
        
        if ~flag_noROspike
            output(triali,6) = round(brady_count);
            
            output(triali,8) = round(preibi,2);
            output(triali,9) = round(post0,2);
            output(triali,10) = round(post1,2);
            output(triali,11) = round(post2,2);
            output(triali,12) = round(post3,2);
            
            output(triali,14) = round(perc0,2);
            output(triali,15) = round(perc1,2);
            output(triali,16) = round(perc2,2);
            output(triali,17) = round(perc3,2);
            output(triali,18) = find(cond_codes == VSDI.condition(triali,1));
        end
        
        
        output(triali,20) = measure(roi2plotidx(1));
        output(triali,21) = measure(roi2plotidx(2));
        
        %     ----------------------------------------------
        %     COPY RESULTS INT A CELL TO PRINT
        %     ------------------------------------------------
        % ----------------------------------------------
        % RESULS - triali-referenced so they can be subset to group
        % ----------------------------------------------
        outputcell{ii,1} = output(triali,1);
        outputcell{ii,2} = output(triali,2);
        outputcell{ii,3} = output(triali,3);
        outputcell{ii,4} = output(triali,4);
        
        if ~flag_noROspike
            outputcell{ii,6} = output(triali,6);
            
            outputcell{ii,8} = output(triali,8);
            outputcell{ii,9} = output(triali,9);
            outputcell{ii,10} = output(triali,10);
            outputcell{ii,11} = output(triali,11);
            outputcell{ii,12} = output(triali,12);
            
            outputcell{ii,14} = output(triali,14);
            outputcell{ii,15} = output(triali,15);
            outputcell{ii,16} = output(triali,16);
            outputcell{ii,17} = output(triali,17);
            outputcell{ii,18} = output(triali,18);

        end
        
        outputcell{ii,20} = round(output(triali,20),3);
        outputcell{ii,21} = round(output(triali,21),3);
        ii = ii+1;
        
        %             clearvars -except VSDI spike movies F0 roi2plotidx output j export w outfolder window method window_sh nfish excelname
        
        
    end %triali
    
    %----------------------------------------------------------------
    % @DEFINE PARAMETERS
    %----------------------------------------------------------------
    paramet{1,1} = 'preStim (s)';
    
    %     paramet{1,1}= 'mA ';
    paramet{1,2}= w.pre;
    paramet{1,3}= 'postStim (s)';
    paramet{1,4}= w.post;
    
    paramet{2,1}= 'wmean';
    
    paramet{3,1} = selroinames{1};
    paramet{3,2} = selroinames{2};

    paramet{4,1}= 'wmean window:'; paramet{4,2}= [num2str( window.wmean(1)) 'to' num2str( window.wmean(2))  'ms'];
    paramet{5,1}= 'baseline:'; paramet{5,2}=  [num2str(window.basel(1)) 'to' num2str(window.basel(2))  'ms'];
    
    paramet{7,1} = date;
    paramet{8,1} = 'source' ;
    paramet{8,5} =mfilename('fullpath');


    % END OF SETTINGS
    
    %     -----------------------------------------------
    %     ADD LABELS IN THE FIRST ROW
    %     ------------------------------------------------
            % write output (new sheet for each fish
        labels = {'idx','triali','cond','spike(s)','','brady_count','','preIbi','post0','post1','post2','post3','','%ibi0' , '%ibi1','%ibi2','%ibi3','','',selroinames{1},selroinames{2}, 'ncond'};  % labels
        for col= 1:numel(labels)
            outputcell{1,col} = labels{col};
        end

    %     -----------------------------------------------
    %     WRITE EXCEL GROUPED
    %     ------------------------------------------------
    
    if export.excel_indiv == 1
        
        %         output2 = cell2table(output2);
       writecell(outputcell, excelname_i, 'sheet', [num2str(VSDI.ref), '_', num2str(cond_codes(end))])

        writecell (paramet, excelname_i, 'sheet', ['param' num2str(VSDI.ref)])
        
    end
    
    if export.excel_singlefish
        writecell(outputcell, excelname_1)
        paramet{9,5} = 'bins from the first fish';
        writecell (paramet, excelname_1, 'sheet', 'params')
    end

    
    % ----------------------------------------------
    % MEAN VALUES FOR EACH CONDITION
    % ----------------------------------------------
    %     triali_kinds = unique(VSDI.condition(:,1));
    %     triali_kinds = triali_kinds(~isnan(triali_kinds));
    
    % GET ONLY SELECTED VALUES
    
    % HEADER
    localoutput{1,1} = 'cod';
    localoutput{1,2} = 'n_cond';
    localoutput{1,3} = 'mA';
    localoutput{1,4} = '%brady (beats)'; %brady count
    localoutput{1,5} = ''; %brady count
    localoutput{1,6} = '%ibi0'; % ibis0 %
    localoutput{1,7} = '%ibi1'; % ibis1 %
    localoutput{1,8} = '%ibi2'; % ibis2 %
    localoutput{1,9} = '%ibi3'; % ibis3 %
    
    jj = 11;
    for roi_i =roi2plotidx
        localoutput{1,jj} = roilabels{roi_i};
        jj = jj+1;
    end
    
    k = 2;
    c = 0;
    for which_cond =  makeRow(cond_codes)
        c = c+1; % n of condition (to get equal condition labels for R)
        
        localoutput{k,1} = which_cond;
        localoutput{k,2} = c;
            mAidx = find(VSDI.condition(:,1) == which_cond, 1,'first');
            mA = VSDI.condition(mAidx,4);
        localoutput{k,3} = mA;%ncond
            trialscond = find(output(:,3)==which_cond); 
        localoutput{k,4} = round(mean(output(trialscond,6)),2); %brady count
        
        localoutput{k,6} = round(mean(output(trialscond,14)),2); % ibis0 %
        localoutput{k,7} = round(mean(output(trialscond,15)),2); % ibis1 %
        localoutput{k,8} = round(mean(output(trialscond,16)),2); % ibis2 %
        localoutput{k,9} = round(mean(output(trialscond,17)),2); % ibis3 %
        
        localoutput{k,11} = round(mean(output(trialscond,20)),3); %
        localoutput{k,12} = round(mean(output(trialscond,21)),3); %
        k = k+1;
        
    end
    
    
    %     figure
    %     plot(1:length(condlabel),roiactiv(:, 1));
    %     hold on
    %     plot(roiactiv(:, 2));
    %     plot(roiactiv(:, 3));
    %     plot(roiactiv(:, 4));
    %     legend dm4 dm4m dm2 dldm
    %
    %     xticklabels(condlabel)
    %
    %     title ([num2str(VSDI.ref), ': peak minus baseline'])
    
    %     ----------------------------------------------
    %     WRITE EXCEL GROUPED
    %     ----------------------------------------------
    if export.excel == 1
        % write output (new sheet for each fish
        %         localoutput = cell2table(localoutput);
        writecell (localoutput, excelname_g, 'sheet', [num2str(VSDI.ref),'_', num2str(cond_codes(end))])
        writecell (paramet, excelname_g, 'sheet', ['param' num2str(VSDI.ref)])
    end

    if export.excel_singlefish
        writecell(localoutput, excelname_1, 'sheet', 'bycond')
    end

    % return
    
    %% ----------------------------------------------
    %     COMPOUND PLOT : (1) CLOUDPLOT + REGRESSION LINE; (2) BINNED SCATTER;
    %     (3) MEANS FROM BINS
    
    %----------------------------------------------------------------
    % @SET: RECALCULATE REJECTION SETTINGS
    %----------------------------------------------------------------
    
    % IBI = 1; % 0 , 1
    % whattrials = 'sharks'; % 'sharks' 'no-sharks' 'all'
    
    % GET SELECTED CONDITIONS TRIALS
    conditrials = [];
    for condi = cond_codes
        conditrials = [conditrials; find(VSDI.condition(:,1) == condi)];
    end
    conditrials = sort(conditrials);
    
    
    % (1) CLOUDPLOT
    % -----------------------------------------------------------
    if IBI == 0
        bradyoutput = output(:, 14);
    elseif IBI==1
        bradyoutput = output(:, 15);
    end
    
    % --------------------------------------------------------------
    % GET TRIALS THAT ARE INSIDE THE SELECTED BRADYCARDIA RANGE
    % --------------------------------------------------------------
    if bradyrange
        
        inrangeidx = [];
        
        for ii = 1:length(bradyoutput)
            if bradyoutput(ii)>=bradyrange(1) && bradyoutput(ii)<=bradyrange(end)
                inrangeidx = [inrangeidx; ii];
            end
        end
        seltrials = intersect(seltrials, inrangeidx);
    end

    
    % --------------------------------------------------------------
    % BINNED BRADY
    % --------------------------------------------------------------
if IBI == 0
brady_all = output(:, 14);
elseif IBI==1
brady_all = output(:, 15);
end 

%build bin-vector
binning_vector = [lim1:step:lim2];
% get binning categories
[~, binId] = histc(brady_all, binning_vector) ;    % Bin edges: 0 --bin1-- 4 --bin2-- 10

% OUT OF RANGE VALUES ARE BINNED INTO '0'

% build labels
for bin = 2:numel(binning_vector)
bin_labels{bin-1} = [num2str(binning_vector(bin-1)) 'to' num2str(binning_vector(bin))];
end

line_roiname = selroinames{1}; % ONLY THE FIRST ROI WILL BE PLOTTED
act1= makeRow(output(:,20)); ...CORRESPONDING TO THIS COL %20 = selroiname{1}; 21=selroiname{2} 

% Get mean activity for each bin
i = 1;
for bini = 1:max(binId)
idx = find(binId == bini);
idx = intersect(idx, seltrials);
meanbin(i) = bini;
meanact1(i) = mean(act1(idx));
i =i+1;
end 

    if export.excel_singlefish
        temp = makeCol(meanbin);
        temp(:,2) = makeCol(meanact1);
        %write separate file
        writematrix(temp, excelname_1, 'sheet', 'bins')
        writecell(bin_labels, excelname_1, 'sheet', 'binlabels')

        clear temp
        
    end


if forceline
   idxnoNan = find(~isnan(meanact1)); 
   if numel(idxnoNan) > 3
       meanact1 = meanact1(idxnoNan);
       meanbin = meanbin(idxnoNan);
   end
end

% ------------------------------------
% Calculating the regression coefficient
% ------------------------------------

% Compute the mean of your data
x= makeRow(meanbin); 

y1 = makeCol(meanact1);

% ref data 
save2plot(rowi).fish = VSDI.ref;
save2plot(rowi).cond = cond_codes;
save2plot(rowi).kind = cond_description;

% data to plot binned means + regression line
save2plot(rowi).bin = x; 
save2plot(rowi).meanact = y1; 
save2plot(rowi).roi = line_roiname; 
save2plot(rowi).binlabels = bin_labels;
 
clearvars -except groupbrady  whattrials  outfield dataunits  save2plot group_note selrow
end


% PLOT FROM VALUES SAVED IN 'save2plot'
figure
for rowi =selrow
    fishref = save2plot(rowi).fish; 
    nfish = TORus('who', fishref);
    VSDI = TORus('load',nfish); 
% ------------------------------------
% Calculating the regression coefficient
% ------------------------------------

x = save2plot(rowi).bin;
y1= save2plot(rowi).meanact;
line_roiname= save2plot(rowi).roi; 
bin_labels= save2plot(rowi).binlabels; 
cond_description= save2plot(rowi).kind;
% ------------------------------------
% Calculating the regression coefficient
% ------------------------------------
% Compute the mean of your data
y1hat = mean(y1);
xhat = mean(x);

% Compute regression coefficients in the least-square sense
b1 = (x - xhat)*(y1 - y1hat)/sum((x - xhat).^2); % Regression coefficient
a1 = y1hat - b1*xhat;                             % Y-intercept

ccmap = lines(7);
ccmap(8,:) = [1 .5 .5]; 
ccmap(9,:)= [.8 .8 .8]; [rho, p] = corrcoef(x, y1);
rho= round(rho(2),1);
p = round(p(2),3);

% ------------------------------------
% Plot results
% ------------------------------------

hold on
scatter(x, y1, 45,  ccmap(rowi,:), 'filled', 'HandleVisibility','off' ) %'Displayname', [num2str(VSDI.ref) '-' cond_description]
plot(x, a1 + b1*x, 'color', ccmap(rowi,:), 'linewidth', 2, 'Displayname', [num2str(VSDI.ref) '-' cond_description  '-' line_roiname '(' num2str(p) ')']); %['rho' num2str(rho) '(p' num2str(p) ')']

clear x y1 line_labels cond_description b1 a1 rho p
end

legend( 'location', 'southwest')
set(gca, 'Xlim', [0 numel(bin_labels)])
set(gca,'xtick',[1:numel(bin_labels)],'xticklabel',bin_labels)
xtickangle(45)
sgtitle(['REGRESSION LINE - ' whattrials])
legend('location', 'southwestoutside')

set(gca,'FontSize',15) % Creates an axes and sets its FontSize to 18
ylabel([line_roiname ':' outfield '(' dataunits ')']);
xlabel(['brady (%ibi)']);
set(gcf, 'Position', get(0, 'Screensize'));

%% Update history:
% 13/04/22 - Created from devo_bradi_ibis_grouped_by_condition_new3_subset_allfish.m
% 06/04/22 - add regression line for all fish; checked correspondence;  fixed bug: idx = intersect(idx, seltrials)
% among: (1) excel of individual trials, (2) excel grouped by conditions,
% (3) points from regression line 
% Created: 05/04/2022 (from: 'devo_bradi_ibis_grouped_by_condition_new2_subset_bradyrange'  Updated on: 31/03/22)