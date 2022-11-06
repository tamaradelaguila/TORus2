% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% OUTPUT TO BE PROCESSED WITH: 
% ~/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/Zscore2/individual_barplot_panel2.R

close all
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus';
user_settings
cd(W)

% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new2.mat') % ATTENTION: select cases manually (there are repetitions)
savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/Zscore2' ;%@ SET

% ///////////////////////////////////////////////////////////
% SETTINGS

selroinames = {'dm4m_R2','dm2_R2'}; %4roi

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_18filt6';% '_17filt5' ; '_18filt6'
% ref_movie= '_12filt5' ;

activ_unit = 'diffF'; % @ MANUALLY SET (just for display purposes)
analysisref = 'intra'; % MANUALLY SET!!! extra info for the name. Set group according to the rows of 'groupplot' selected in;   for suji  =  [1 3 4 9] lateral_group7_RECOV
reject_on = 3;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/Zscore2' ;%@ SET

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 1200]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-25 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350];

% feedf.noise.fr_abovenoise = 30;
% feedf.noise.SDfactor = 2;
% feedf.noise.basel = [-200 0]; %it won't be used by the function, but will be used to manually get the baseline
feedf.method = 'movsum';

% flagslope = 1;

% Params
slope.window = [0, 200]; %ms

% END USER_SETTINGS
% ///////////////////////////////////////////////////////////////

% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{4,1} = 'wmean (ms)';
params{4,2} = feedf.window.wmean;

params{5,1} = 'slope win (ms)';
params{5,2} = slope.window;


if strcmpi( analysisref, 'panel2_group8')
    sel_subjects = [1 2 5 6];
end

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

% Subsettings:
setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET+
setting.force_include = 0; %@ SET


%% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new2.mat') % ATTENTION: select cases manually (there are repetitions)

% suji = 5;
% nfish = groupplot{suji,1};
% cond_list = groupplot{suji,3};

nfish = 11;
cond_list = [400 401 402 404];

VSDI = TORus('load', nfish);
VSDmov = TORus('loadmovie',nfish,ref_movie);
movies = VSDmov.data ;
F0 = VSDmov.F0;

%----------------------------------------------------------------
% GET INDEXES OF TIMERANGE
%----------------------------------------------------------------
idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values

timebase_adj = VSDI.timebase(idxrange);

% get indexes for slope
slope.windowidx = dsearchn(timebase_adj, makeCol(slope.window));
slope.windowidx = [slope.windowidx(1) slope.windowidx(end)];

%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------
rej = 'reject' ;
if reject_on > 1
    rej = [rej num2str(reject_on)];
end

rejectidx = [];

if setting.manual_reject
    try
        rejectidx = [rejectidx  makeRow(VSDI.(rej).manual)];
    catch
        rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
        disp(['reject.manual was used for fish' num2str(VSDI.ref) 'because there is no reject' num2str(reject_on) '.manual'])
    end
end

if setting.GSabsthres_reject
    rejectidx = [rejectidx  makeRow(VSDI.(rej).GSabs025)];
    
end

if setting.GSmethod_reject
    rejectidx = [rejectidx makeRow(VSDI.(rej).GSdeviat2sd)];
end

if setting.force_include
    rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
    
end

rejectidx = sort(unique(rejectidx));


%         params{8,1} = 'rejected';

%         params{8,1} = rejectidx';

%----------------------------------------------------------------
% SELECT ROI
%----------------------------------------------------------------
switch roikind
    case 'circle'
        selroi =name2idx(selroinames, VSDI.roi.labels_circ);
        roilabels = VSDI.roi.labels_circ;
        masks =  VSDI.roi.circle.mask;
        
    case 'anat'
        selroi =name2idx(selroinames, VSDI.roi.labels);
        roilabels = VSDI.roi.labels;
        masks = VSDI.roi.manual_mask;
end

% -------------------------------------------
% GET INDIVIDUAL MEASURES
% -------------------------------------------
rowi= 1 ;

for condition =  makeRow(cond_list)
    % -------------------------------------------
    % SELECT CASES
    % -------------------------------------------
    sel_trials = find(VSDI.condition(:,1)==condition);
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    
    for triali = makeRow(sel_trials)
        trialmov = movies(:,:,idxrange,triali);
        
        for roi_i = makeRow(selroi)
            roimask= masks(:,:,roi_i);
            trialwave= roi_TSave(trialmov,roimask);
            
            % ---------------------------------------------------------------
            ... CALCULATE MEASURES FROM TRIAL
            % ---------------------------------------------------------------
            
            temp = devo_peak2peak(trialwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
            
            roipeak = temp.peakminusbasel;
            roiwmean= temp.wmean;
            
            %slope
            idx0 = slope.windowidx(1);
            idxend = slope.windowidx(end);
            waveW = trialwave(idx0:idxend);
            slopemean = mean(diff(waveW));
            roislope = slopemean;
            clear pixelwave waveW slopemean
            
            % ---------------------------------------------------------------
            ... STORE INTO LONG-FORMAT CELL
                % ---------------------------------------------------------------
            rowi = rowi+1; %first row has to be 2 (first is for labels)
            % FACTORS
            longF{rowi,1} = VSDI.ref ; %subject id
            longF{rowi,2} = roi_i; % roi
            longF{rowi,3} = roilabels{roi_i}; % roi
            longF{rowi,4} = condition; % condition (1 = low; 2 = hi; 3 = blank)
            longF{rowi,5} = triali; % condition (1 = low; 2 = hi; 3 = blank)
            % MEASURES
            longF{rowi,6} = round(roipeak,2); % outputP(roi_i, condi)
            longF{rowi,7} = round(roiwmean,2);%
            longF{rowi,8} = round(roislope,2);%
            
        end
    end % for triali
    
end %for condition


%% -------------------------------------------
% EXPORT FOR R
% -------------------------------------------
excelname = fullfile(savein, [num2str(VSDI.ref) 'long_format_forR' ref_movie '_rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls']);
labels = {'fish' 'roin' 'roi' 'cond' 'trialidx' 'peak' 'wmean' 'slope'};

% Assign labels to first row
for col = 1:numel(labels)
    longF{1,col}= labels{col};
end

% write output (new sheet for each fish
writecell (longF, excelname, 'sheet',  [roikind ref_movie 'rej' num2str(reject_on)])

writecell (params, excelname, 'sheet','param')

clear longF

blob()

% Created: 24/02/22 from srouce code:
% /home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/03_figure_sketch/groupplot_measures/group_plot_brain_reprogram_4_Zscore2_working.m