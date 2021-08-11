%% test1_settings

% user_settings

   path1=  '/home/tamara/Documents/MATLAB/VSDI/TORus';
   addpath(genpath(path1))
   addpath(genpath('/home/tamara/Documents/MATLAB/VSDI/VSDI_ourToolbox'))
   
   path.data = fullfile(path1, 'data');
    path.grouplist = fullfile(path.rootpath);
    path.list =fullfile(path.rootpath, 'data','BVlists');

    
    %% TEST PARAMETERS SETTINGS

    wind.min = [-100 100];
    wind.max = [0 600];
    wind.movsum = 50;
    
    nroi = 3; 
   
    pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_peak2peak_function';
    
    method_kind = {'movsum', 'findpeaks'};
    cond_codes = [2000 2001 2002 2003]; %for trial-wise wave

   