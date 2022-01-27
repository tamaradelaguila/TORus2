%% SEED-BASED connectivity
% fast analysis 


%% MANUAL SEED DEFINITION
% from the 


clear

%----------------------------------------------------------------
% @SET: BASIC PARAMETERS
%----------------------------------------------------------------

ref_movie= '_12filt5' ;
% ref_movie= '_11filt4' ;
% ref_movie = '_08diff_perc_f0pre';
% ref_movie = '_09filt3';

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')


%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

reject_on = [12];  %@ SET
% Subsettings:
setting.manual_reject = 0; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET+
setting.force_include = 0; %@ SET


% END OF SETTINGS


% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

% end of user_settings

%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------


 blocki = 8;% 1:length(fast_condition_list)
   
    nfish = fast_condition_list{blocki,1};
    
    trial_kinds = fast_condition_list{blocki,2};
    cond_def =fast_condition_list{blocki,3};
    
    
    [VSDI] = TORus('load',nfish);
    
    temp = TORus('loadmovie',nfish,ref_movie);
    movies = temp.data(:,:,1:end-1,:);

%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------

rejectidx = [];

if setting.manual_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
end

if setting.GSabsthres_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];

end

if setting.GSmethod_reject
    rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];
end

if setting.force_include
    rejectidx = setdiff(rejectidx, VSDI.reject.forcein);

end

rejectidx = sort(unique(rejectidx));

%----------------------------------------------------------------
% SELECT CASES
%----------------------------------------------------------------
sel_trials= [];

for condi = trial_kinds %to make sure that only the conditions of interest are computed)
    condtrials = makeCol(find(VSDI.condition(:,1)==condi));
    sel_trials  = [sel_trials; condtrials];
end
sel_trials= sort(sel_trials);

if reject_on  %@ SET
    sel_trials = setdiff(sel_trials, rejectidx);
end



%
%----------------------------------------------------------------
% DRAW ROI (first time) OR INPUT COORDINATES 
%----------------------------------------------------------------
  
 [cROI] = roicirc_draw1ext(VSDI, 3);
    

%% OPTION 1 - CONCATENATE MOVIES
case 'concat_simple'
%----------------------------------------------------------------
% CONCATENATE TRIALS FROM EACH CONDITION 
%----------------------------------------------------------------
j = 1;

for condi =  makeRow(trial_kinds)
   idxcondi =  find(VSDI.condition(:,1)==condi);
   idxcondi = intersect(idxcondi, sel_trials);
%    idxcondi = intersect(idxcondi(1:round(end/2)), sel_trials);
%    idxcondi = intersect(idxcondi(round(end/2):end), sel_trials);
   
moviesin = movies(:,:,:,idxcondi);
dim = size(moviesin);
newdim = dim(3)*dim(4);

supermovie = reshape(moviesin, [dim(1) dim(2) newdim]);

seedwave = roi_TSave(supermovie,cROI.mask);
seedwave = seedwave';

%----------------------------------------------------------------
% PIXEL-WISE CORRELATION 
%----------------------------------------------------------------
for rowi =1:dim(2)
    for coli = 1:dim(1)
        pixelwave = squeeze(supermovie(coli, rowi, 1:end-1)); %roi_TSave substracts the last frame so we have to match length
        
        if sum(pixelwave) ==0
            rhomap(coli, rowi, j) = 0;
        else
            rhomap(coli, rowi, j) = corr(seedwave, pixelwave); 
        end
        
        clear pixelwave
    end
end
j = j+1;
clear idx moviesin
end
% 
% 
% case 'concat_match'
% 
%     j = 1;
% 
% for condi =  makeRow(trial_kinds)
%    idxcondi =  find(VSDI.condition(:,1)==condi);
%    idxcondi = intersect(idxcondi, sel_trials);
% %    idxcondi = intersect(idxcondi(1:round(end/2)), sel_trials);
% %    idxcondi = intersect(idxcondi(round(end/2):end), sel_trials);
%    
% moviesin = movies(:,:,:,idxcondi);
% dim = size(moviesin);
% 
% for triali = 1:dim(4)
% supermovie = reshape(moviesin, [dim(1) dim(2) newdim]);
% end 
% 
% seedwave = roi_TSave(supermovie,cROI.mask);
% seedwave = seedwave';
% 
% %-------
% % pixel-wise correlation 
% %----------------------------------------------------------------
% for rowi =1:dim(2)
%     for coli = 1:dim(1)
%         pixelwave = squeeze(supermovie(coli, rowi, 1:end-1)); %roi_TSave substracts the last frame so we have to match length
%         
%         if sum(pixelwave) ==0
%             rhomap(coli, rowi, j) = 0;
%         else
%             rhomap(coli, rowi, j) = corr(seedwave, pixelwave); 
%         end
%         
%         clear pixelwave
%     end
% end
% j = j+1;
% clear idx moviesin
% end

case 'trialwise'

    
end % cases method

%----------------------------------------------------------------
% PLOT 
%----------------------------------------------------------------

% GET mA corresponding to each condition to plot in the graph
j =1 ;
for condi =  makeRow(trial_kinds)
idxcondi =  find(VSDI.condition(:,1)==condi, 1);
mA(j) = VSDI.condition(idxcondi,4);
j = j+1;
end

% PLOT THE CORRELATION (RHO) MAPS
figure
nplot = length(trial_kinds)+1;
for ploti = 1:nplot-1
    subplot(1,nplot,ploti);
    imagesc(rhomap(:,:,ploti))
    colormap(flipud(hot))
    
    axis image
    set(gca, 'clim', [0 1])
    title(['cond=' num2str(mA(ploti)) 'mA'])
%     freezeColors()
     colorbar
%     cbfreeze(h(ploti))

end

% PLOT SEED ROI IN THE LAST PLOT

ax1= subplot(1,nplot,nplot);
imagesc(VSDI.crop.preview); hold on
colormap(ax1, 'bone')
axis image
title('seed')
viscircles(cROI.centre,cROI.r, 'color','r')
axis image

supertit = [num2str(VSDI.ref)  'seed-based correlation'];
sgtitle(supertit)

name = [VSDI.ref '_data_' ref_movie 'seed-based correlation'];

%% OPTION 2 - CONCATENATE MOVIES WITHOUT THE BASELINE
%----------------------------------------------------------------
% CONCATENATE TRIALS FROM EACH CONDITION 
%----------------------------------------------------------------
j = 1;
analysislabel = 'concatenated';

for condi =  makeRow(trial_kinds)
   idxcondi =  find(VSDI.condition(:,1)==condi);
%    idxcondi = intersect(idxcondi, sel_trials);
%    idxcondi = intersect(idxcondi(1:round(end/2)), sel_trials);
   idxcondi = intersect(idxcondi(round(end/2):end), sel_trials);
   
   
   basetimeidx = dsearchn(VSDI.timebase, 0) ;
   
moviesin = movies(:,:,basetimeidx:end,idxcondi);
dim = size(moviesin);
newdim = dim(3)*dim(4);

supermovie = reshape(moviesin, [dim(1) dim(2) newdim]);

seedwave = roi_TSave(supermovie,cROI.mask);
seedwave = seedwave';

%----------------------------------------------------------------
% PIXEL-WISE CORRELATION 
%----------------------------------------------------------------
for rowi =1:dim(2)
    for coli = 1:dim(1)
        pixelwave = squeeze(supermovie(coli, rowi, 1:end-1)); %roi_TSave substracts the last frame so we have to match length
        
        if sum(pixelwave) ==0
            rhomap(coli, rowi, j) = 0;
        else
            rhomap(coli, rowi, j) = corr(seedwave, pixelwave); 
        end
        
        clear pixelwave
    end
end
j = j+1;
clear idx moviesin
end


%----------------------------------------------------------------
% PLOT 
%----------------------------------------------------------------

% GET mA corresponding to each condition to plot in the graph
j =1 ;
for condi =  makeRow(trial_kinds)
idxcondi =  find(VSDI.condition(:,1)==condi, 1);
mA(j) = VSDI.condition(idxcondi,4);
j = j+1;
end

% PLOT THE CORRELATION (RHO) MAPS
figure
nplot = length(trial_kinds)+1;
for ploti = 1:nplot-1
    subplot(1,nplot,ploti);
    imagesc(rhomap(:,:,ploti))
    colormap (flipud(hot))
    

    axis image
    set(gca, 'clim', [0 1])
    title(['cond=' num2str(mA(ploti)) 'mA'])
%     freezeColors()
     colorbar
%     cbfreeze(h(ploti))

end

% PLOT SEED ROI IN THE LAST PLOT

ax1= subplot(1,nplot,nplot);
imagesc(VSDI.crop.preview); hold on
colormap(ax1, 'bone')
axis image
title('seed')
viscircles(cROI.centre,cROI.r, 'color','r')
axis image

supertit = [num2str(VSDI.ref)  'seed-based correlation (no baseline)'];
sgtitle(supertit)

name = [VSDI.ref '_data_' ref_movie 'seed-based correlation no baseline'];

%% OPTION 3 - CORRELATION TRIALWISE AND THEN AVERAGED
j = 1;

for condi =  makeRow(trial_kinds)
   idxcondi =  find(VSDI.condition(:,1)==condi);
%    idxcondi = intersect(idxcondi, sel_trials);
%    idxcondi = intersect(idxcondi(1:round(end/2)), sel_trials);
   idxcondi = intersect(idxcondi(round(end/2):end), sel_trials);
   
   
   basetimeidx = dsearchn(VSDI.timebase, 0) ;
   
moviesin = movies(:,:,basetimeidx:end,idxcondi);
dim = size(moviesin);
newdim = dim(3)*dim(4);

supermovie = reshape(moviesin, [dim(1) dim(2) newdim]);

seedwave = roi_TSave(supermovie,cROI.mask);
seedwave = seedwave';

%----------------------------------------------------------------
% PIXEL-WISE CORRELATION 
%----------------------------------------------------------------
for rowi =1:dim(2)
    for coli = 1:dim(1)
        pixelwave = squeeze(supermovie(coli, rowi, 1:end-1)); %roi_TSave substracts the last frame so we have to match length
        
        if sum(pixelwave) ==0
            rhomap(coli, rowi, j) = 0;
        else
            rhomap(coli, rowi, j) = corr(seedwave, pixelwave); 
        end
        
        clear pixelwave
    end
end
j = j+1;
clear idx moviesin
end


%----------------------------------------------------------------
% PLOT 
%----------------------------------------------------------------

% GET mA corresponding to each condition to plot in the graph
j =1 ;
for condi =  makeRow(trial_kinds)
idxcondi =  find(VSDI.condition(:,1)==condi, 1);
mA(j) = VSDI.condition(idxcondi,4);
j = j+1;
end

% PLOT THE CORRELATION (RHO) MAPS
figure
nplot = length(trial_kinds)+1;
for ploti = 1:nplot-1
    subplot(1,nplot,ploti);
    imagesc(rhomap(:,:,ploti))
    colormap (flipud(hot))
    

    axis image
    set(gca, 'clim', [0 1])
    title(['cond=' num2str(mA(ploti)) 'mA'])
%     freezeColors()
     colorbar
%     cbfreeze(h(ploti))

end

% PLOT SEED ROI IN THE LAST PLOT

ax1= subplot(1,nplot,nplot);
imagesc(VSDI.crop.preview); hold on
colormap(ax1, 'bone')
axis image
title('seed')
viscircles(cROI.centre,cROI.r, 'color','r')
axis image

supertit = [num2str(VSDI.ref)  'seed-based correlation (from trialw)'];
sgtitle(supertit)

name = [VSDI.ref '_data_' ref_movie 'seed-based correlation trialw ave'];


%% PLOT CONDITION MINUS BASELINE


%----------------------------------------------------------------
% PLOT 
%----------------------------------------------------------------

% GET mA corresponding to each condition to plot in the graph
j =1 ;
for condi =  makeRow(trial_kinds)
idxcondi =  find(VSDI.condition(:,1)==condi, 1);
mA(j) = VSDI.condition(idxcondi,4);
j = j+1;
end

% PLOT THE CORRELATION (RHO) MAPS
figure
nplot = length(trial_kinds);
for ploti = 1:nplot-1
    subplot(1,nplot,ploti);
    cond = ploti+1; %we skip the baseline map
    imagesc(rhomap(:,:,cond) - rhomap(:,:,1)) % substract the first condition
    colormap (polarmap())
    

    axis image
    set(gca, 'clim', [-0.5 0.5])
    title(['cond=' num2str(mA(cond)) 'mA'])
%     freezeColors()
     colorbar
%     cbfreeze(h(ploti))

end

% PLOT SEED ROI IN THE LAST PLOT

ax1= subplot(1,nplot,nplot);
imagesc(VSDI.crop.preview); hold on
colormap(ax1, 'bone')
axis image
title('seed')
viscircles(cROI.centre,cROI.r, 'color','r')
axis image

supertit = [num2str(VSDI.ref)  'differential seed-based correlation (cond minus control)'];
sgtitle(supertit)

name = [VSDI.ref '_data_' ref_movie 'differential seed-based correlation'];