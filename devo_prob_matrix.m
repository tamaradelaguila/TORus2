%% 
clear
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')

% ///////////////////////////////////////////////////////////
% SETTINGS

selroinames = {'dm4m_R', 'dm4_R',  'dm2_R', 'dm3_R', 'dldm_R', 'dm4m_L', 'dm4_L', 'dm2_L', 'dm3_L', 'dldm_L',};
% roikind = 'circle'; %
roikind = 'anat';

ref_movie= '_17filt5' ;
% ref_movie= '_12filt5' ;
% ///////////////////////////////////////////////////////////


%  'SETTINGS'
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');

addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));
% END 'USER_SETTINGS'

%% GET TRANSITION MATRICES 
ref_wave = 'circ_filt618';

flagcirc = strcmpi(ref_wave(1:4), 'circ');
    

%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
for  block = 1%[1:4]% 1:length(fast_condition_list)
    
    % Get selection to be analyzed from structure
%     nfish = groupplot{block,1};
nfish = 11;
%     trial_kinds = groupplot{block,3};
trial_kinds = [400 404] ;
    
    [VSDI] = TORus('load',nfish);
    
    if flagcirc
    selroi = name2idx(selroinames, VSDI.roi.labels_circ);
    else 
    selroi = name2idx(selroinames, VSDI.roi.labels);
    end
    
        selroi = [roi1 roi2];

    
    VSDroiTS = TORus('loadwave',nfish);
    waves = VSDroiTS.(ref_wave).data; %@ SET
    
end 

%% PLOT INTO BIOGRAPH