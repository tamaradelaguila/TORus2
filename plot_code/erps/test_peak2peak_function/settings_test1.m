
% user settings

path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; 
 tempsep.idcs   = strfind(path.rootpath,'/');
 tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

clear temp

%% TEST SETTINGS


    wind.min = [-100 100];
    wind.max = [0 600];
    wind.movsum = 50;
    
    nroi = 3; 
   
    pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_peak2peak_function';
    
    method_kind = {'movsum', 'findpeaks'};
    cond_codes = [2000 2001 2002 2003]; %for trial-wise wave
