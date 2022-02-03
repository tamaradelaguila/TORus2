% CODE by Marco Marino 

clear all
close all
clc

load('data_example_tamara.mat')

% name variables
filtered = data.filtered;
unfiltered = data.unfiltered;

% x*y*time*trials
size(filtered)
size(unfiltered)

% check data
% figure;
% subplot(211)
% imagesc(squeeze(mean(unfiltered(:,:,200,:),4)))
% caxis([-5 5])
% subplot(212)
% imagesc(squeeze(mean(filtered(:,:,200,:),4)))
% caxis([-5 5])

data_unfiltered = mean(unfiltered(:,:,:,:),4);
data_filtered = mean(filtered(:,:,:,:),4);

% reshape data before applying PCA or ICA
for i = 1:681
    
data_tmp = squeeze(data_unfiltered(:,:,i));
data_tmp = data_tmp(:);
data_tmp = data_tmp';
data_tot(i,:) = data_tmp;

end

% data_tot = detrend(data_tot);

randn('seed',0);
[COEFF, SCORE, LATENT] = pca(data_tot,'Centered','off');
LATENT=100*LATENT/sum(LATENT); % Percentage of explained variance associated with each component

% reconstruct all PCs spatial pattern
for i = 1:10 % for the first 10 components for example

data_tmp = COEFF(:,i);
B = reshape(data_tmp,[60 89]);
pattern_pcs(:,:,i) = B;

figure(i);
subplot(211)
imagesc(squeeze(pattern_pcs(:,:,i)))
subplot(212)
plot(SCORE(:,i))

end

% you can select how many PCs to include (this can be implemented also in an automated way)
n_pcs = [1:10] % e.g. from the first to the fifth component (excluding the second)
data_pcs_tmp = SCORE(:,n_pcs)*COEFF(:,n_pcs)';

% back reconstruct the images (matrix) in time from the vector for each
% time point
for i = 1:681

data_tmp = data_pcs_tmp(i,:);
B = reshape(data_tmp,[60 89]);
data_pcs(:,:,i) = B;

end

t = 200 % time frame
figure;
subplot(311)
imagesc(squeeze(data_unfiltered(:,:,t)))
caxis([-8 8]) % specify the range of the colorbar
subplot(312)
imagesc(squeeze(data_filtered(:,:,t)))
caxis([-8 8]) % specify the range of the colorbar
subplot(313)
imagesc(squeeze(data_pcs(:,:,t)))
caxis([-8 8]) % specify the range of the colorbar

coord_x = 30;
coord_y = 20;

figure;
plot(squeeze(data_unfiltered(coord_x,coord_y,:)),'k')
hold on
plot(squeeze(data_filtered(coord_x,coord_y,:)),'b')
hold on
plot(squeeze(data_pcs(coord_x,coord_y,:)),'r')

%% INDEPENDENT COMPONENT ANALYSIS
n_comp = 10; % the number of components expected should be selected

% spatial ICA time (t) x pixels (p) [t x p]
% ICA on unfiltered data
% [IC, mixing, unmixing]=fastica(data_tot,'approach','defl','numOfIC',n_comp,'g','tanh','maxNumIterations',500);

addpath(genpath('/home/tamara/Documents/MATLAB/VSDI/TORus/ICA'))
% ICA on unfiltered data where PCA was applied
[IC, mixing, unmixing]=fastica(data_pcs_tmp,'approach','defl','numOfIC',n_comp,'g','tanh','maxNumIterations',500);
% IC=unmixing*data_tot;

% back reconstruct the images (matrix) in time from the vector for each
% time point
for i = 1:n_comp

data_tmp = IC(i,:);
B = reshape(data_tmp,[60 89]);
data_IC(:,:,i) = B;

figure(20+i);
subplot(211)
imagesc(B)
caxis([-5 5])
subplot(212)
plot(mixing(:,i))

end

% select the good components from the figures you have just plotted (look at the spatial pattern and at the time-course)
good_ics = [1,3,6,7]

data_ics_tmp = mixing(:,good_ics)*IC(good_ics,:);

for i = 1:681

data_tmp = data_ics_tmp(i,:);
B = reshape(data_tmp,[60 89]);
data_IC(:,:,i) = B;

end

t = 200
figure;
subplot(311)
imagesc(squeeze(data_unfiltered(:,:,t)))
caxis([-8 8]) % specify the range of the colorbar
subplot(312)
imagesc(squeeze(data_filtered(:,:,t)))
caxis([-8 8]) % specify the range of the colorbar
subplot(313)
imagesc(squeeze(data_IC(:,:,t)))
caxis([-8 8]) % specify the range of the colorbar


% MEAN MAPS 
figure;
t1= 100; t2 = 300;
subplot(311)
imagesc(squeeze(mean(data_unfiltered(:,:,t1:t2),3)))
% caxis([-8 8]) % specify the range of the colorbar
subplot(312)
imagesc(squeeze(mean(data_filtered(:,:,t1:t2),3)))
% caxis([-8 8]) % specify the range of the colorbar
subplot(313)
imagesc(squeeze(mean(data_IC(:,:,t1:t2),3)))
% caxis([-8 8]) % specify the range of the colorbar
% 

coord_x = 20;
coord_y = 30;

figure;
plot(squeeze(data_unfiltered(coord_x,coord_y,:)),'k')
hold on
plot(squeeze(data_filtered(coord_x,coord_y,:)),'b')
hold on
plot(squeeze(data_IC(coord_x,coord_y,:)),'r')

% 
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
% data_tmp = mixing(:,i);
% B = reshape(data_tmp,[60 89]);
% data_IC(:,:,i) = B;
% 
% figure(i);
% imagesc(B)
% % caxis([-0.1 0.1])
% end

