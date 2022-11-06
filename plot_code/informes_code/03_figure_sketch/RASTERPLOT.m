clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus';
user_settings
cd(W)

%----------------------------------------------------------------
% @SET: BASIC PARAMETERS
%----------------------------------------------------------------

% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new.mat') % ATTENTION: select cases manually (there are repetitions)

%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
nfish = 8;
cond_codes = [1000];

ref_movie= '_18filt6' ;

saveraster = 0;
rasterclim = [] ; %when [], the script is useful for outliers identification

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/rasterplot'; %@ SET

selroinames = {'dm4m_R2',  'dm2_R2' ,'dm1_R','dldm_R2','dm4m_L2', 'dm2_L2' ,'dm1_L','dldm_L2'};%dm3
roikind = 'circle'; %

trange = [-300 1434]; %ms Range of analysis

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on= 3;
setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.force_include = 0; %

%% LOAD / COMPUTE SETTINGS
%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus ('load', nfish);
VSDmov = TORus('loadmovie',nfish,ref_movie);

%----------------------------------------------------------------
... GET INDEXES OF TIMERANGE AND ADJUSTED TIMEBASE
%----------------------------------------------------------------
idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values

idx0 = dsearchn(makeCol(VSDI.timebase), 0);
timebase_adj = VSDI.timebase(idxrange);

%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------
rej = 'reject' ;
if reject_on > 1
    rej = [rej num2str(reject_on)];
end

rejectidx = [];

if setting.manual_reject
    rejectidx = [rejectidx  makeRow(VSDI.(rej).manual)];
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

        %----------------------------------------------------------------
        % ADJUST SELECTED ROI  (according to whether the roi exists in the
        % fish
        %----------------------------------------------------------------
%         ii = 1;
%         selroinames = [];
%         for roii = 1:length(selroinames_ini)
%             idx =  find(contains( VSDI.roi.labels,selroinames_ini{roii}));
%             if sum(idx)>0
%                 selroinames{ii} = VSDI.roi.labels{idx};
%                 ii = ii+1;
%             end
%             clear idx
%         end
% 


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


%% ----------------------------------------------------------------
... BUILD RASTER
    %----------------------------------------------------------------
for condi = makeRow(cond_codes)
[sel_trials] = find(VSDI.condition(:,1)==condi);

if reject_on
    sel_trials = setdiff(sel_trials, rejectidx);
    disp('trials rejected')
end

n = numel(timebase_adj)-1;
roiraster = NaN(numel(sel_trials), n,  numel(selroi));

ii = 0; %counter for all rasters in 2D

for roii = makeRow(selroi)
    
    roimask = masks(:,:,roii);
    
    tri = 0;
    for triali = makeRow( sel_trials)
        tri = tri +1;
        ii = ii+1;
        %to plot single trial
        movie2plot = squeeze(VSDmov.data(:,:,idxrange,triali));
        meanF0 = squeeze(VSDmov.F0(:,:,triali));
        
        roiraster(tri,:,roii) =  roi_TSave_percF_roiwise(movie2plot,roimask, meanF0);
        rasters2D (ii,1:n) =  roi_TSave_percF_roiwise(movie2plot,roimask, meanF0);
    end % for triali
    ii = ii +1; 
    rasters2D (ii:ii+1,:) =  NaN(2,n);
    ii = ii +2; 
end % for roii

%% ----------------------------------------------------------------
... PLOT RASTER
%------------------------------------------------------------------
maxval = max(roiraster(:));
rasterclim = [0 maxval*0.9];
nroi = numel(selroi);

figure
ploti = 0;
for roii = makeRow(selroi)
    ploti = ploti+1;
    h(ploti) = subplot(nroi,1,ploti);
    imagesc(roiraster(:,:,roii));
    axis tight
    if isempty(rasterclim)
        set(h(ploti),'xtick',[],'ytick',[])
        disp('ATT: raster color limits adjusted to its own maximum')
    else
        set(h(ploti),'xtick',[],'ytick',[], 'clim', rasterclim)
    end
    
    xline(idx0, 'color', 'w');
    ylabel(roilabels{roii})
    colormap(jet)
    
    savename= ['RoiRasters_' num2str(VSDI.ref) ref_movie '_cond' num2str(condi)  '_rej'  num2str(reject_on) '.jpg'];
    titulo = [num2str(VSDI.ref) '.All trials from cond:' num2str(condi) ];
    
    if saveraster
        %           saveas(gcf, fullfile(savein, savename), 'jpg')
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
        print(fullfile(savein, [savename]),'-r600','-djpeg') % prints it as you see them
        close
    end
    
    
end

end % for condi
