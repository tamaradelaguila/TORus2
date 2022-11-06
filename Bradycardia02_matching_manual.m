%% BRADYCHARDIA MATCHING CODE - MANUALLY MATCHED PREVIOUSLY

% Bradycardia has to be prevously calculated from the spike file with the
% counting code BRADYCARDIA_flexible_counting_code.m). The previous steps needed are defined in  'spike_export.txt'

% In this code, both the matching trial and the bradychardia will be
% extracted from the same input CSV, that will be obtained as in the following description: 

% Manually copy the columns from the output excel of 'Bradycardia01_flexible_counting' into a CSV that will be the input for the code. 
% The CSV has to contain, therefore: 
% (1) A 'stim' column 
% (2) A 'trial' column with the BrainVision reference's number of the trial (eg: 003). 
% Other columns (with variable name) containing the brady measures for the
% trials

% So, 'stim' and 'trials' need to have fixed names, unlike
% The rest of columns with the brady measures, whose names will be retrieve
% in the code to be stored as they were named

%% 
clear
nfish = 13; 
user_settings

% Csv with brady data (*data manually extracted from the output excell from the counting code)
sourcebrady.in = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike/excels_brady'; 
sourcebrady.filelist= '220611spike_BRADYcodeinput.csv';

VSDI = TORus('load',nfish);


%% STEP 3 (OPTION 1 - FROM CSV ) - MATCHING PROCESS ITSELF
% ....................................................................
% Both, the matching trial and the bradychardia will be extracted from the
% same CSV.
%The CSV can be extracted manually copying columns from the output excel of 'BRADYCARDIA1_flexible_counting' and has to contain: 
% (1) A 'stim' column 
% (2) A 'trial' column with the number of the trial (eg: 003). 
% Other columns (with variable name) containing the brady measures for the
% trials

% So, 'stim' and 'trials' need to have fixed names, unlike
% the rest of columns with the brady measures, whose names will be retrieve in the code


filelistpath = fullfile(sourcebrady.in,sourcebrady.filelist);
% refpath = fullfile( source.in, [num2str(fishref) 'spikeref.csv']);

% IMPORT CVS FILE
bradydata = readtable(filelistpath,'ReadVariableNames',true, 'delimiter', ',');
bradydata = table2struct(bradydata);

% TRANSFORM TRIAL NAME TO MATCH THE EXACT BV FILENAME AND FIND THE
% CORRESPONDING IDX !!! run only once 
temp = (VSDI.ref * 1000); 
for ii = 1: length(bradydata)
    triali =  bradydata(ii).trial + temp;
    bradydata(ii).trial = triali; 
    bradydata(ii).idx = find(VSDI.trialref == triali);
    
    if isempty(bradydata(ii).idx ) % it's important to fill the [] with NaN to avoid errors when using 'find([bradydata(:).idx] == value)'
    bradydata(ii).idx = NaN;
    end
end

% % UNCOMPLETE - ALTERNATIVE WITHOUT LOOP
% temp = [tablelist2.trial] + (VSDI.ref * 1000); %get vector from the structure and perform computation
% % BUT... how to assign the new vector straight to the structure???? 


% the loop has to advance respect to the structure, because if we extract
% the vector from the structure, empty values will be skipped and therefore
% the idx changed

% GET LIST OF BRADY MEASURES CONTAINED ('stim' and 'trial' are remove and
% the rest are brady measures)
bradylist = fieldnames(bradydata);
bradylist = setdiff(bradylist, {'stim'}); 
bradylist = setdiff(bradylist, {'trial'}); 
bradylist = setdiff(bradylist, {'idx'}); 

% SAVE BOTH STRUCTURES
VSDI.brady.data = bradydata;
VSDI.brady.measurelist = bradylist;

% TORus('save', VSDI)

%% TO FURTHER USE THE BRADY MEASURES

% To find a specific trial measure: first find the idx of the bradydata
idx = find([VSDI.brady.data(:).idx] == 33);

% loop through bradylist
for ii= 1:length(VSDI.brady.measurelist) 
field =  bradylist{ii};
bradi = VSDI.brady.data.(field); 
end


%% UPDATE HISTORY
% 01/09/22 - Separated code for 'manual' and 'automatic' matching 
% 24/08/22 Finish code for Matching option 2 - manual 
% 19/08/22  Created from: spike_processing_easyfeedback
