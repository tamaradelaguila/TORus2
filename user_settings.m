% if 'pwd' is used instead of manually adding the rootpath, user_settings
% has to be executed everytime from the rootpath folder (make sure that you
% are in the rootpath folder for your experiment everytime it's going to be
% called from inside a script)



% path.rootpath = 'C:\Users\User\Documents\UGent_brugge\VSDI_tamaraToolbox\VSDI_ourToolbox\rootpath';
path.rootpath = pwd; 
 tempsep.idcs   = strfind(path.rootpath,'/');
 tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');

addpath(genpath(path.rootpath));
addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

clear temp
