% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% BLANK SUBSTRACTION OF THE WAVES

% average movie
% extract wave for each pixel and make the measure maps (for all
% conditions): for condition for xi for yi

% z-spatial all conditions for each fish
% extract roi from z-spatially maps >>> definite measures
close all
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus';
user_settings
cd(W)


% ///////////////////////////////////////////////////////////
% SETTINGS
nfish =13;
cond_codes =[0 301 302];

% selroinames = {'dm4m_L',  'dm2_L'};
selroinames = {'dm4m_R2',  'dm2_R2'};

reject_on = 0;

% selroinames = {'dm4m_R2','dm2_R2', 'dm1_R','dldm_R'}; %4roi

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_18filt6';% 
activ_unit = 'dF'; % @ MANUALLY SET (just for display purposes)


get_excel_indiv = 1; %not working at the moment

savein= '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady2/input';

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
slope.window = [0, 1200]; %ms

% END USER_SETTINGS
% ///////////////////////////////////////////////////////////////

%% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
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

% params{6,1} = ['onset lat ' num2str(onset_factor*100) '%(ms)'];

params{8,1} = date;
params{8,3} = ['source:'  mfilename('fullpath')];


%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

% Subsettings:
setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 0; %@ SET+
setting.forcein = 0; %@ SET


%% CODE: GET INDIVIDUAL  MEASURES OF ROI WAVES FOR EACH TRIAL
    
    VSDI = TORus('load', nfish);
    VSDmov = TORus('loadmovie',nfish,ref_movie);
    movies = VSDmov.data ;
    F0 = VSDmov.F0;
    
    %----------------------------------------------------------------
    % CONTROL ROI PICTURE
    %----------------------------------------------------------------
    %     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R);
    %     title([num2str(VSDI.ref) 'roi preview:' selroinames{:}])
    %     saveas(gcf, [num2str(VSDI.ref)'roipreview'], 'jpg')
    
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
% ROI SELECTION
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

    
%----------------------------------------------------------------
% TRIALS SELECTION - to avoid computing unnecesary trials
%----------------------------------------------------------------
                sel_trials = [];
                for condition =  makeRow(cond_codes)
                    newtri = makeRow(find(VSDI.condition(:,1)==condition));
                    sel_trials = [sel_trials newtri];
                    clear newtri
                end
                
                if reject_on
                    sel_trials = setdiff(sel_trials, rejectidx);
                end

                
%----------------------------------------------------------------
% TRIAL-WISE COMPUTATION 
%----------------------------------------------------------------
                    
                ti = 1;
                
                    for triali = makeRow(sel_trials) %each trial will be saved in a row
                        ti = ti+1;
                        movie2plot = movies(:,:,idxrange,triali);
                        ri= 4;%roi counter
                        
                        for roii = makeRow(selroi) %each roi in a column
                            ri= ri+1;
                            roimask = masks(:,:, roii);
                            
                            %                             meanF0 = squeeze(mean(VSDmov.F0(:,:,triali),3));
                            roiwave = roi_TSave(movie2plot,roimask);
                            fs = 1000/(VSDI.info.stime);
                            roiwave = lowpass(roiwave,5,fs);
                            temp = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window, [], feedf.method, 0, 0);
    
                            ntrial = VSDI.trialref(triali) - VSDI.ref*1000; 
                            
                            % trial identification / kind
                            trialmeasure.peak{ti,1} = triali;
                            trialmeasure.peak{ti,2} = ntrial;
                            trialmeasure.peak{ti,3} = VSDI.condition(triali,1);
                            trialmeasure.peak{ti,4} = VSDI.condition(triali,4);
                            % measures: PEAK
                            trialmeasure.peak{ti,ri} = round(temp.peakminusbasel,2);
    
                            % trial identification / kind
                            trialmeasure.wmean{ti,1} = triali;
                            trialmeasure.wmean{ti,2} = ntrial;
                            trialmeasure.wmean{ti,3} = VSDI.condition(triali,1);
                            trialmeasure.wmean{ti,4} = VSDI.condition(triali,4);
                            % measures: WMEAN
                            trialmeasure.wmean{ti,ri} = round(temp.wmean,2);
    
                            clear roimask roiwave meanF0
    
%                             % TRIALS TABLE TO LATER REJECT BASED ON STD
%                             % trial identification / kind
%                             trialmat.peak(ti,1) = triali;
%                             trialmat.peak(ti,2) = VSDI.trialref(triali);
%                             trialmat.peak(ti,3) = condition;
%                             % measures: PEAK
%                             trialmat.peak(ti,ri) = round(temp.peakminusbasel,2);

                        end % forroii
                    end % for triali
                    
                    % Build labels 
                                    % GET WAVES
                labels = {'idx' 'trial' 'cond' 'mA'};
                ri = 4;
                for roii = makeRow(selroi) %each roi in a column
                    ri = ri+1;
                    labels{1,ri} = roilabels{roii};
                end
                    
                % add labels
                trialmeasure.peak(1,1:length(labels)) = labels(1:end);
                trialmeasure.wmean(1,1:length(labels)) = labels(1:end);
    
                % EXPORT TO 'R'
                excelname = ['corr02_ACTIVITY_fromBV_' num2str(VSDI.ref) ref_movie '_rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls'];
    
                % write output (new sheet for each fish
                writecell (trialmeasure.peak, fullfile(savein , excelname), 'sheet', 'peak')
                writecell (trialmeasure.wmean, fullfile(savein,  excelname), 'sheet', 'wmean')
    
                writecell (params, fullfile(savein,  excelname), 'sheet','param')

                clear trialmeasure
    
   
% Update history:
% 14/11/22: Created (adapted from
% group_plot_brain_reprogram_4_Zscore_working.m)
% 16/02/22 - Add onset-latency measure and delete 'switch refcase' (so only
% the normal dF condition is considered). To use blank-substraction
% methods, go to the older versions of the code
