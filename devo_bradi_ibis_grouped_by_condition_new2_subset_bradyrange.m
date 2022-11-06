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

% SET FISH AND CONDITION TO STUDY
nfish = [1];
cond_codes = [300:303];

% CHOOSE MOVIE REFERENCE ACCORDING TO THE DATAUNITS
dataunits = '%dF'; % '%dF' 'dF'
switch dataunits
    case 'dF'
        ref_movie = '_18filt6'; % input movie '_17filt5'
    case '%dF'
        ref_movie = '_21filt6'; 
end

% EXPORTING SETTINGS
export.excel = 0; % '1' = sí; '0'=no
export.excel_indiv = 0;

%SETTINGS FOR PART [1] OF THE CODE ..............................
w.pre= 10; % ventana pre-estímulo(en segundos) , de los que se va a tomar las espigas para hacer la media de frecuencia instantánea o el conteo
w.post = w.pre/2; % sólo para el conteo

outfield = 'wmean'; % 'peakminusbasel' 'wmean'

outfolder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady';
excelname_i= fullfile(outfolder,[outfield ref_movie '_brady_grouped_bycondition_inrange.xlsx']);
excelname_g =  fullfile(outfolder,[outfield ref_movie '_brady_trialiwise_inrange.xlsx']);

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
selroinames = {'dm4m_R2','dm2_R2'}; % ONLY 2 ROI
% selroinames = {'dm4m_L','dm2_L'}; % just for fish 3I

% ................................................................

% SETTINGS FOR PART [2] OF THE CODE: TRIALS SELECTION FOR PLOTS (do not affect the excel)...
IBI = 1; % 0 , 1
whattrials = 'all'; % 'all', 'sharks', 'no-sharks' 
bradyrange = [-10 60]; %bradychardia range that will be plot (and used to do the regression line)
% if bradyrange = [], the whole range will be used
saveplots = 0; 

% for binned cloudplot
step = 10;
lim1 = -10; 
lim2 = 60;

% ....................................................................

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
seltrials = [];
for condi = cond_codes
    ci = ci+1;
    
    seltrials = [seltrials; find(VSDI.condition(:,1) == condi)];
end

seltrials = sort(seltrials);

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

paramet{1,5}= 'wmean';

paramet{6,1}= 'wmean window:'; paramet{6,2}= [num2str( window.wmean(1)) 'to' num2str( window.wmean(2))  'ms'];
paramet{7,1}= 'baseline:'; paramet{7,2}=  [num2str(window.basel(1)) 'to' num2str(window.basel(2))  'ms'];


paramet{8,4} = 'time window for shark finding spans the whole triali (unlike for roi peak activity)';

paramet{10,1} = selroinames{1};
paramet{10,2} = selroinames{2};
% END OF SETTINGS


%     -----------------------------------------------
%     WRITE EXCEL GROUPED
%     ------------------------------------------------

if export.excel_indiv == 1
    
    % write output (new sheet for each fish
    labels = {'idx','triali','cond','spike(s)','','brady_count','','preIbi','post0','post1','post2','post3','','%ibi0' , '%ibi1','%ibi2','%ibi3','','',selroinames{1},selroinames{2}};  % labels
    for col= 1:numel(labels)
        outputcell{1,col} = labels{col};
    end
    %         output2 = cell2table(output2);
    writecell (paramet, excelname_g, 'sheet', ['param' num2str(VSDI.ref)])
    
    writecell(outputcell, excelname_g, 'sheet', num2str(VSDI.ref))
    
end


% ----------------------------------------------
% MEAN VALUES FOR EACH CONDITION
% ----------------------------------------------
%     triali_kinds = unique(VSDI.condition(:,1));
%     triali_kinds = triali_kinds(~isnan(triali_kinds));

% HEADER
localoutput{1,1} = 'cod';
localoutput{1,2} = 'mA';

localoutput{1,3} = '%brady (beats)'; %brady count

localoutput{1,5} = '%ibi0'; % ibis0 %
localoutput{1,6} = '%ibi1'; % ibis1 %
localoutput{1,7} = '%ibi2'; % ibis2 %
localoutput{1,8} = '%ibi3'; % ibis3 %

jj = 9;
for roi_i =roi2plotidx
    jj = jj+1;
    localoutput{1,jj} = roilabels{roi_i};
end

k = 2;
for which_cond =  makeRow(cond_codes)
    
    localoutput{k,1} = which_cond;
    
    cond_trialis = find(VSDI.condition(:,1) == which_cond);
    
    localoutput{k,2} = VSDI.condition(cond_trialis(1),4);
    
    localoutput{k,3} = round(mean(output(cond_trialis,6)),2); %brady count
    
    localoutput{k,5} = round(mean(output(cond_trialis,14)),2); % ibis0 %
    localoutput{k,6} = round(mean(output(cond_trialis,15)),2); % ibis1 %
    localoutput{k,7} = round(mean(output(cond_trialis,16)),2); % ibis2 %
    localoutput{k,8} = round(mean(output(cond_trialis,17)),2); % ibis3 %
    
    localoutput{k,10} = round(mean(output(cond_trialis,20)),3); %
    localoutput{k,11} = round(mean(output(cond_trialis,21)),3); %
    
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
    writecell (localoutput, excelname, 'sheet', num2str(VSDI.ref))
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

 
        
        
% ADJUST TO THE TRIALS-INCLUSION CRITERIA 
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


% (1) CLOUDPLOT + REGRESSION LINE
% -----------------------------------------------------------
if IBI == 0
bradyoutput = output(:, 14);
elseif IBI==1
bradyoutput = output(:, 15);
end 


% GET TRIALS THAT ARE INSIDE THE SELECTED BRADYCARDIA RANGE

if bradyrange
    
    inrangeidx = [];
    
    for ii = 1:length(bradyoutput)
        if bradyoutput(ii)>=bradyrange(1) && bradyoutput(ii)<=bradyrange(end)
            inrangeidx = [inrangeidx; ii];
        end
    end
    seltrials = intersect(seltrials, inrangeidx);
end


brady = bradyoutput(seltrials); 
act1= output(seltrials,20)';
act2= output(seltrials,21)';

% Calculating the regression coefficient
% Compute the mean of your data
x= brady'; 

y1hat = mean(act1);
y2hat = mean(act2);

xhat = mean(x);

% Compute regression coefficients in the least-square sense
b1 = (x - xhat)*(act1' - y1hat)/sum((x - xhat).^2); % Regression coefficient
b2 = (x - xhat)*(act2' - y2hat)/sum((x - xhat).^2); % Regression coefficient

a1 = y1hat - b1*xhat;                             % Y-intercept
a2 = y2hat - b2*xhat;                             % Y-intercept

figure
subplot(1,3,1)
scatter(x, act1, 40,  'r', 'filled')
hold on
scatter(brady, act2, 40,  'b', 'filled')

plot(x, a1 + b1*x, 'r', 'linewidth', 1);
plot(x, a2 + b2*x, 'b', 'linewidth', 1);

legend([selroinames(1) , selroinames(2)], 'location', 'northwest')
ylabel([outfield '(' dataunits ')']);
xlabel(['brady (%ibi' num2str(IBI) ')']);
title(['cloudplot + regression line'])

set(gca,'FontSize',15) % Creates an axes and sets its FontSize to 18


% plotregression(brady, act1, 'Regression')

% (2)+(3) BINNED MEASURE
% -----------------------------------------------------------
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

% Get bins and activity from selected truals 
binId_sel = binId(seltrials);
act1= output(seltrials,20)';
act2= output(seltrials,21)';

% Get mean activity for each bin
i = 1;
for bini = min(binId):max(binId)
idx = find(binId_sel == bini);
% idx = intersect(idx, seltrials);

meanbin(i) = bini;
meanact1(i) = mean(act1(idx));
meanact2(i) = mean(act2(idx));
i =i+1;
end 


% Plot
subplot(1,3,2) % binned scatter 
scatter(binId_sel, act1, 40, 'r')
hold on
scatter(binId_sel, act2, 40,  'b')
scatter(meanbin, meanact1, 70, 'd', 'r', 'filled')
scatter(meanbin, meanact2, 70, 'd', 'b', 'filled')

legend([selroinames(1), selroinames(2)], 'location', 'northwest')
set(gca,'xtick',[1:numel(bin_labels)],'xticklabel',bin_labels)
ylabel([outfield '(' dataunits ')']);
xtickangle(45)
newX = get(gca, 'XLim');
% set(gca, 'Xlim', [newX(1)-1 newX(2)+1]) 
set(gca, 'Xlim', [newX(1) newX(2)+1])  % aavoid plotting '0' (out of range values)
xlabel(['brady (%ibi' num2str(IBI) ')']);
title(['binned cloudplot + mean'])

set(gca,'FontSize',15) % Creates an axes and sets its FontSize to 18


subplot(1,3,3) % means tendency
scatter(meanbin, meanact1, 70, 'd', 'r', 'filled')
hold on
scatter(meanbin, meanact2, 70, 'd', 'b', 'filled')

legend([selroinames(1), selroinames(2)], 'location', 'northwest')
set(gca,'xtick',[1:numel(bin_labels)],'xticklabel',bin_labels)
ylabel([outfield '(' dataunits ')']);
xtickangle(45)
newX = get(gca, 'XLim');
set(gca, 'Xlim', [newX(1) newX(2)+1])
xlabel(['brady (%ibi' num2str(IBI) ')']);
title('means from bins')

sgtitle(['fish' num2str(VSDI.ref) '- (cond' num2str(cond_codes) ')' '--' whattrials])

set(gca,'FontSize',15) % Creates an axes and sets its FontSize to 18

clear brady act1 act2 x xhat y1hat y2hat a1 a2 b1 b2 x step lim1 lim2 binning_vector binId bin_labels binId_sel meanbin
% return
pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady'; 
side = selroinames{1}(end);
name1 = ['SHARKvsnoSHARK_INRANGE' num2str(VSDI.ref), '_', num2str(cond_codes(1)), 'ibi', num2str(IBI), '_', whattrials, '_reject', num2str(reject_on),'_', side, '_BRADI.jpg' ];
set(gcf, 'Position', get(0, 'Screensize'));
if saveplots
    saveas(gcf, fullfile(pathsave, name1), 'jpg')
    close
end
%% DRAW WAVES 

% waveroinames = {'dm4m_R2','dm2_R2', 'dm4m_L2', 'dm2_L', 'dldm_R'}; 
waveroinames = {'dm4m_R','dm2_R', 'dm4m_L', 'dm2_L', 'dldm_R'}; 

switch roikind
    case 'circle'
        roi2wave =name2idx(waveroinames, VSDI.roi.labels_circ);
        roiwavelabels = VSDI.roi.labels_circ;
        masks =  VSDI.roi.circle.mask;
        
    case 'anat'
        roi2wave =name2idx(waveroinames, VSDI.roi.labels);
        roiwavelabels = VSDI.roi.labels;
        masks = VSDI.roi.manual_mask;
end

%----------------------------------------------------------------
% SELECT triali
%----------------------------------------------------------------
ncond = length(cond_codes);

figure
sgtitle(['fish' num2str(VSDI.ref)  '. Trials:' whattrials])

% PLOT ROI PREVIEW

        cmap = roicolors();
        cmap = cmap(1:2:end,:); %roicolors map has double values for 2 hemispheres
        
        back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                
        sp1= subplot(2,ncond,1)
                
        switch roikind
            case 'circle'
                centers = VSDI.roi.circle.center(roi2wave, :) ;
                roicirc_preview_multiple_cmap_HIdef(back, centers, VSDI.roi.circle.R, sp1, cmap,5);
                

            case 'anat'
                roi_preview_multiple(back, VSDI.roi.manual_poly(roi2wave,:), sp1);
        end

            
cond_codes_no0 = [];
for condi = cond_codes
strcode = num2str(condi);
if strcmpi(strcode(end),'0')
    cond_codes_no0= cond_codes_no0;
else
    cond_codes_no0 = [cond_codes_no0 condi];
end
   
end
    
ci = 0;
for condi = cond_codes_no0
    ci = ci+1;
    
    subplot(2,ncond,ci+ncond)
    title(num2str(condi))
    % ----------------------------------------------
    % GET WAVES AND PLOT
    % ----------------------------------------------
    conditrials = find(VSDI.condition(:,1) == condi); 
    trials = intersect(conditrials, seltrials); 
    
    movie = movies(:,:,:,trials);
    movie = mean(movie,4); 
    
    for nroi = makeRow(roi2wave)
        roimask= masks(:,:,nroi);
        
        roiwave =  roi_TSave(movie,roimask);
        plot(VSDI.timebase, roiwave, 'linewidth', 1.2)
        hold on
    end    
    clear movie conditrials trials
end
set(gcf, 'Position', get(0, 'Screensize'));

if saveplots
    name2 = ['SHARKvsnoSHARK_INRANGE' num2str(VSDI.ref), '_', num2str(cond_codes(1)), 'ibi', num2str(IBI), '_', whattrials, '_reject', num2str(reject_on),'_', side, '.jpg' ];
    saveas(gcf, fullfile(pathsave, name2), 'jpg')
    close
end


%% ----------------------------------------------
%     ALL FISH: CLOUDPLOT : (1) CLOUDPLOT + REGRESSION LINE; 

%----------------------------------------------------------------
% @SET: RECALCULATE REJECTION SETTINGS
%----------------------------------------------------------------

load()
% IBI = 1; % 0 , 1
% whattrials = 'sharks'; % 'sharks' 'no-sharks' 'all'

% GET SELECTED CONDITIONS TRIALS 
        conditrials = [];
        for condi = cond_codes
            conditrials = [conditrials; find(VSDI.condition(:,1) == condi)];
        end
        conditrials = sort(conditrials);

        
% ADJUST TO THE TRIALS-INCLUSION CRITERIA 
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


% (1) CLOUDPLOT + REGRESSION LINE
% -----------------------------------------------------------
if IBI == 0
bradyoutput = output(:, 14);
elseif IBI==1
bradyoutput = output(:, 15);
end 


% GET TRIALS THAT ARE INSIDE THE SELECTED BRADYCARDIA RANGE

if bradyrange
    
    inrangeidx = [];
    
    for ii = 1:length(bradyoutput)
        if bradyoutput(ii)>=bradyrange(1) && bradyoutput(ii)<=bradyrange(end)
            inrangeidx = [inrangeidx; ii];
        end
    end
    seltrials = intersect(seltrials, inrangeidx);
end


brady = bradyoutput(seltrials); 
act1= output(seltrials,20)';
act2= output(seltrials,21)';

% Calculating the regression coefficient
% Compute the mean of your data
x= brady'; 

y1hat = mean(act1);
y2hat = mean(act2);

xhat = mean(x);

% Compute regression coefficients in the least-square sense
b1 = (x - xhat)*(act1' - y1hat)/sum((x - xhat).^2); % Regression coefficient
b2 = (x - xhat)*(act2' - y2hat)/sum((x - xhat).^2); % Regression coefficient

a1 = y1hat - b1*xhat;                             % Y-intercept
a2 = y2hat - b2*xhat;                             % Y-intercept

figure
subplot(1,3,1)
scatter(x, act1, 40,  'r', 'filled')
hold on
scatter(brady, act2, 40,  'b', 'filled')

plot(x, a1 + b1*x, 'r', 'linewidth', 1);
plot(x, a2 + b2*x, 'b', 'linewidth', 1);

legend([selroinames(1) , selroinames(2)], 'location', 'northwest')
ylabel([outfield '(' dataunits ')']);
xlabel(['brady (%ibi' num2str(IBI) ')']);
title(['cloudplot + regression line'])

set(gca,'FontSize',15) % Creates an axes and sets its FontSize to 18

%% Created: 14/07/2021
% Update history:
% 06/04/22: fixed bug: idx = intersect(idx, seltrials);
% Add that bradyrange = [] computes the whole range; add pearson coef
% 31/03/22: add trials subselection (for part 2 of the code): 'all',
% 'sharks', 'no-sharks'
% 30/03/22: add code to draw waves
% 29/03/22: refine 'noisewin set
% 24/03/22: add dataunits. Record 'ref_movie' parameter in the name
% 08/03/22 add binning of data to plot in ranges
