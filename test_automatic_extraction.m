%% CHECK IMPORT .MAT FILE EXPORTED FROM PYTHON ARRAY
% EXPLORE NUMPY.ARRAY EXPORTED WITH THE PYTHON COMMANDS: 

% Import from python
load('C:\Users\User\Documents\UGent_brugge\lab_maps\BVdml\output\200611000A.mat')
dataPy= double(data_act); dataPy = permute(dataPy, [2 3 1]); %to match dimensions from BVimport
dataPy = dataPy(1:59,1:87,:); %THE FIRST IMAGE IS THE BACKGROUND! 
% EXTRACTED VALUES: 
% 1 - frame0 (background)
% 2:681 - frames 1 a 680 (6a4080ms)

% Import from BV 
[originalBV,~]= BVimport('C:\Users\User\Documents\UGent_brugge\lab_maps\BVdml\200611000A-waveEach.csv', 680);
dataBV(:,:,1) = originalBV(:,:,682);%leave background out
dataBV(:,:,2:681) = originalBV(:,:,1:680);

minval = min(dataPy(:));
dataTam = dataPy + abs(minval);

% Extract vector
x=20; y=25; trial = 1;
vectorPy = squeeze(dataPy(x,y,:));
vectorBV = squeeze(dataBV(x,y,:));

vectorTam =  squeeze(dataTam(x,y,:));

plot(vectorPy); hold on; plot(vectorBV); plot(vectorTam); legend py BV tam

corr(vectorBV, vectorTam)


