% CODE by Marco Marino

clear all
close all
clc

user_settings

nfish = 11;
cond = 403;
reject_on =0; %4 (4 for connectivity analysis)

% ----------------------------------------------------------------
... @SET: REJECT PARAMETERS
% ----------------------------------------------------------------

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.force_include = 0; %

% ----------------------------------------------
...  LOAD DATA
%     ----------------------------------------------

VSDI = TORus('load', nfish);
VSDmov = TORus('loadmovie', nfish, '_16diff_f0pre');
VSDmov_filt =  TORus('loadmovie', nfish, '_18filt6');

% ----------------------------------------------------------------
... COMPUTE REJECTION IDX
% ----------------------------------------------------------------
rej = 'reject' ;
if reject_on > 1
    rej = [rej num2str(reject_on)];
end

rejectidx = [];


if setting.manual_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
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

% ----------------------------------------------------------------
... IDX SELECTION
% ----------------------------------------------------------------

idxsel = choose_condidx(VSDI.condition(:,1),cond);

if reject_on
    idxsel = setdiff(idxsel, rejectidx);
end

idxsel = sort(idxsel');

% ----------------------------------------------------------------
... SELECT DATA 
% ----------------------------------------------------------------

filtered = VSDmov_filt.data(:,:,1:end-1,idxsel) ; 
unfiltered = VSDmov.data(:,:,1:end-1,idxsel) ; 



% Get dimensions
d= size(filtered);
nrow = d(1);
ncol = d(2);
ntime = d(3);

% x*y*time*trials
% size(filtered)
% size(unfiltered)

% VISUALY CHECK DATA
figure;
subplot(211)
imagesc(squeeze(mean(unfiltered(:,:,200,:),4)))
caxis([-5 5])
subplot(212)
imagesc(squeeze(mean(filtered(:,:,200,:),4)))
caxis([-5 5])

% ----------------------------------------------------------------
... AVERAGE TRIALS 
% ----------------------------------------------------------------

data_unfiltered = mean(unfiltered(:,:,:,:),4);
data_filtered = mean(filtered(:,:,:,:),4);


% ALTERNATIVE: CONCATENATED TRIALS
rawdim = size(unfiltered);
data_unfiltered = reshape(unfiltered, [rawdim(1), rawdim(2),rawdim(3)*rawdim(4)]);
data_filtered = reshape(filtered, [rawdim(1), rawdim(2),rawdim(3)*rawdim(4)]);
d = size(data_unfiltered);


% ----------------------------------------------
...  RESHAPE DATA for PCA and ICA : time * pix
% ----------------------------------------------

% for i = 1:ntime
%     data_tmp = squeeze(data_unfiltered(:,:,i));
%     data_tmp = data_tmp(:);
%     data_tmp = data_tmp';
%     data_tot(i,:) = data_tmp; % pix*timepix*pix*timetime
%     
% end

data_tot = reshape(data_unfiltered, [d(1)*d(2) , d(3)]); 
data_tot =data_tot'; 

% ----------------------------------------------
...  PCA and RECONSTRUCTION
% ----------------------------------------------

% data_tot = detrend(data_tot);

randn('seed',0);
[COEFF, SCORE, LATENT] = pca(data_tot,'Centered','off');
LATENT=100*LATENT/sum(LATENT); % Percentage of explained variance associated with each component

% reconstruct all PCs spatial pattern

for coefi = 1:10 % for the first 10 components for example
    figure
    subplot(212)
    
    data_tmp = COEFF(:,coefi);
    data_reconstr = reshape(data_tmp,[nrow ncol]);
    
    subplot(211)
    imagesc(squeeze(data_reconstr))
    
    subplot(212)
    plot(SCORE(:,coefi))
    
    sgtitle(['PC=' num2str(coefi)])
    pattern_pcs(:,:,coefi) = data_reconstr;
    
end

%% ----------------------------------------------
...  SELECT PCs TO INCLUDE
%     ----------------------------------------------
% you can select how many PCs to include (this can be implemented also in an automated way)
n_pcs = [1:5]; % e.g. from the first to the fifth component (excluding the second)
data_pcs_tmp = SCORE(:,n_pcs)*COEFF(:,n_pcs)';

% ----------------------------------------------
...  MOVIE RECONSTRUCTION
%     ----------------------------------------------
% back reconstruct the images (matrix) in time from the vector for each
% time point
for ti = 1:ntime
    
    data_tmp = data_pcs_tmp(ti,:);
    B = reshape(data_tmp,[nrow ncol]);
    data_pcs(:,:,ti) = B;
    
end

% VISUAL CHECK
t = 200;  % time frame
figure;
subplot(311)
title('raw data')
imagesc(squeeze(data_unfiltered(:,:,t)))
colorbar
% caxis([-8 8]) % specify the range of the colorbar
subplot(312)
title('filt data')
imagesc(squeeze(data_filtered(:,:,t)))
colorbar
% caxis([-8 8]) % specify the range of the colorbar
subplot(313)
imagesc(squeeze(data_pcs(:,:,t)))
title('PC-reconstr from raw')
colorbar
% caxis([-8 8]) % specify the range of the colorbar

coord_x = 30;
coord_y = 20;

figure;
plot(squeeze(data_unfiltered(coord_x,coord_y,:)),'k')
hold on
plot(squeeze(data_filtered(coord_x,coord_y,:)),'b')
hold on
plot(squeeze(data_pcs(coord_x,coord_y,:)),'r')
legend 'raw' 'filt' 'pc'

% INDEPENDENT COMPONENT ANALYSIS
n_comp = 10; % the number of components expected should be selected (set to 10)

% spatial ICA time (t) x pixels (p) [t x p]
% ICA on unfiltered data
% [IC, mixing, unmixing]=fastica(data_tot,'approach','defl','numOfIC',n_comp,'g','tanh','maxNumIterations',500);

% ICA on unfiltered data where PCA was applied
[IC, mixing, unmixing]=fastica(data_pcs_tmp,'approach','defl','numOfIC',n_comp,'g','tanh','maxNumIterations',500);
% IC=unmixing*data_tot;

% BACK RECONSTRUCTION FROM COMPONENTS AND VISUALIZATION
for i = 1:size(IC,1)% n_comp
    
    data_tmp = IC(i,:);
    B = reshape(data_tmp,[nrow ncol]);
    data_IC(:,:,i) = B;
    
    figure;
    subplot(211)
    imagesc(B)
    colorbar
    caxis([-5 5])
    subplot(212)
    plot(mixing(:,i))
    
    sgtitle(['IC=' num2str(i)]) 
end

%% select the good components from the figures you have just plotted (look at the spatial pattern and at the time-course)
good_ics = [1 2 3 5 6];

data_ics_tmp = mixing(:,good_ics)*IC(good_ics,:);

for i = 1:ntime
    
    data_tmp = data_ics_tmp(i,:);
    B = reshape(data_tmp,[nrow ncol]);
    data_IC(:,:,i) = B;
    
end

t = 200;
figure;
subplot(311)
imagesc(squeeze(data_unfiltered(:,:,t)))
colorbar
title('raw data')
% caxis([-8 8]) % specify the range of the colorbar
subplot(312)
imagesc(squeeze(data_filtered(:,:,t)))
colorbar
title('filt data')
% caxis([-8 8]) % specify the range of the colorbar
subplot(313)
imagesc(squeeze(data_IC(:,:,t)))
colorbar
title('IC reconstr (from raw data)')
% caxis([-8 8]) % specify the range of the colorbar

% MEAN MAPS
figure;
t1= 100; t2 = 300;
subplot(311)
imagesc(squeeze(mean(data_unfiltered(:,:,t1:t2),3)))
colorbar
title('raw data')
% caxis([-8 8]) % specify the range of the colorbar
subplot(312)
imagesc(squeeze(mean(data_filtered(:,:,t1:t2),3)))
colorbar
title('filt data')
% caxis([-8 8]) % specify the range of the colorbar
subplot(313)
imagesc(squeeze(mean(data_IC(:,:,t1:t2),3)))
colorbar
title('IC reconstr (from raw data)')
% caxis([-8 8]) % specify the range of the colorbar


coord_x = 20;
coord_y = 30;

figure;
plot(squeeze(data_unfiltered(coord_x,coord_y,:)),'k')
hold on
plot(squeeze(data_filtered(coord_x,coord_y,:)),'b')
hold on
plot(squeeze(data_IC(coord_x,coord_y,:)),'r')
legend 'raw' 'filt' 'IC'

% % temporal ICA pixels (p) x time (t) [p x t]
% [IC, mixing, unmixing]=fastica(data_tot','approach','defl','numOfIC',n_comp,'g','tanh','maxNumIterations',500);
% % IC=unmixing*data_tot;
% 
% % figure;
% % imagesc(mix)
% % data_ics_tmp = mixing(:,good_ics)*IC(good_ics,:);
% 
% % back reconstruct the images (matrix) in time from the vector for each
% % time point
% for i = 1:n_comp
%     
%     data_tmp = mixing(:,i);
%     B = reshape(data_tmp,[60 89]);
%     data_IC(:,:,i) = B;
%     
%     figure(i);
%     imagesc(B)
%     % caxis([-0.1 0.1])
% end
% 
