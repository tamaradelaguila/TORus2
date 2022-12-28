%% DRAW CIRCULAR ROI WITH FUNCTION
clear
user_settings

% movie to preview wave
ref_movie= '_18filt6';% '_15filt5', '_17filt5' ; '_18filt6'

nfish = 31; 
condition = 303;

reject_on =0;

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis


preview_roilabels = {'dm4m_R2',...
    'dm2_R2',...
    'dm1_R',...
    'dldm_R'};


preview_roilabels = {'dm4m_R2', 'dm2_R2'};

%% ----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus('load',nfish);
VSDmov = TORus('loadmovie',nfish,ref_movie);
movies = VSDmov.data ;
F0 = VSDmov.F0; % when present

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

% Subsettings:
setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 0; %@ SET+
setting.forcein = 0; %@ SET

% ----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------

rejectidx = [];
rejectidx  = compute_rejectidx(VSDI, reject_on, setting);


    % -------------------------------------------
    % SELECT CASES  AND AVERAGE MOVIE
    % -------------------------------------------
    sel_trials = find(VSDI.condition(:,1)==condition);

    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end

    %----------------------------------------------------------------
% GET INDEXES OF TIMERANGE
%----------------------------------------------------------------
idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values

timebase_adj = VSDI.timebase(idxrange);

    
    % --------------------------------------------------------------------------
    % GET AVERAGE MOVIE. Use timerange set
    % --------------------------------------------------------------------------

    movieave = mean(movies(:,:,idxrange,sel_trials),4);
    movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
    % (that normally corresponds to the background)

%% ROI LABELS

newlabel{1} = 'new_dm2_R2';
close all

%STEP 1: select R and draw roi

r = 4; % r = 4 for whole brain
[coord, mask] = roicir_draw(VSDI.crop.preview,newlabel,r); 


%% STEP 2 - preview
circle.center = coord;
circle.R = r;
circle.mask = mask;

 roicirc_preview_multiple(VSDI.crop.preview, circle.center, circle.R); 

 %% STEP 3 - CHECK IN WAVE PLOT (TOGETHER WITH OTHER REGIONS)
 
 % WAVES IN 'dF'
 figure
 subplot(1,2,1)
 hold on
 % get preview waves
 for ii = 1:numel(preview_roilabels)
     name_loc = preview_roilabels{ii};
     idx = name2idx(name_loc,VSDI.roi.labels_circ);
     mask = VSDI.roi.circle.mask(:,:,idx);
     wave_local_dF = roi_TSave(movieave,mask);
     wave_local_dF = movmean(wave_local_dF ,5);
     plot(timebase_adj,wave_local_dF, 'linewidth',2, 'displayname',  name_loc); hold on
 end
 
 % get new wave
 roiwave_dF = roi_TSave(movieave,circle.mask);
 roiwave_dF = movmean(roiwave_dF ,5);
 plot(timebase_adj,roiwave_dF, 'linewidth',2, 'displayname',  newlabel{1}); hold on
 legend
 
 hold off
 ax= subplot(1,2,2);
 roicirc_preview_multiple(VSDI.crop.preview, circle.center, circle.R, ax);
 
 sgtitle([num2str(VSDI.ref) '- cond' num2str(condition), '(units diffF).'  newlabel{1} '=' num2str(round(circle.center))])

% WAVES IN 'pF'

F0ave = F0(:,:,sel_trials);
figure
subplot(1,2,1)

% get preview waves
for ii = 1:numel(preview_roilabels)
    name_loc = preview_roilabels{ii};
    idx = name2idx(name_loc,VSDI.roi.labels_circ);
    mask = VSDI.roi.circle.mask(:,:,idx);
    wave_local_pF= roi_TSave_percF_roiwise(movieave,mask, F0ave);
    wave_local_pF = movmean(wave_local_pF ,5);
    plot(timebase_adj,wave_local_pF, 'linewidth',2, 'displayname',  name_loc); hold on
end

% get new wave
roiwave_pF = roi_TSave_percF_roiwise(movieave,circle.mask, F0ave);
plot(timebase_adj,roiwave_pF, 'linewidth',2, 'displayname',  newlabel{1}); hold on
hold off
legend

ax= subplot(1,2,2);
roicirc_preview_multiple(VSDI.crop.preview, circle.center, circle.R, ax);

sgtitle([num2str(VSDI.ref) '- cond' num2str(condition), '(units %F).'  newlabel{1} '=' num2str(round(circle.center))])


    ... if you don't like it, go to step 1
     ... if you like it, save in to STEP 4
 

%% STEP 4: assign a roi number (idx) to store and save:
newidx = 16;

VSDI.roi.labels_circ{newidx} = newlabel{1};
VSDI.roi.circle.center(newidx,:) = circle.center;
% VSDI.roi.circle.R = circle.R;
VSDI.roi.circle.mask(:,:,newidx) = circle.mask;
TORus('save',VSDI)
