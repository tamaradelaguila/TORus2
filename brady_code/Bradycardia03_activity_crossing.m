%% CODE TO CROSS DATA FROM BRADI + ACTIVATION FROM CODES:
% bradi data from spike: 
% activity data from trial's waves: /home/tamara/Documents/MATLAB/VSDI/TORus/brady_code/Individual_trial_measures.m
clear

folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/brady2/input';

file.spike  = 'corr01_BRADY_fromSpike_220615correlacion.xls';
sheet.spike = 'ordered';

file.ROmatlab  = 'corr02_ACTIVITY_fromBV_220615_18filt6_rej0_2roi.xls';
sheet.ROmatlab = 'wmean';

% check that there was no error
ref_spike = file.spike(24:29);
ref_RO = file.ROmatlab(24:29);

if ~strcmpi(ref_spike, ref_RO)
    error('check that the fish is the same')
end

% GET ACTIVITY, CONDITION, TRIAL FROM MATLAB EXCEL

[RO.data, RO.labels] = xlsread(fullfile(folder, file.ROmatlab), sheet.ROmatlab); 

temp = find(strcmpi([RO.labels(:)], {'trial'}));
trial.RO = RO.data(:,temp);
[spike.data, spike.labels] = xlsread(fullfile(folder, file.spike), sheet.spike); 
spike.labels = spike.labels(1,:);

temp = find(strcmpi([spike.labels(:)], {'bvtrial'}));
trial.spike = spike.data(:,temp);
clear temp

% GET TIMES OF SPIKE AND BRADI VECTORS
time = spike.data(:,1);


temp = find(strcmpi([spike.labels(:)], {'%ibi'}));
bradi = spike.data(:,temp);
clear temp

temp = find(strcmpi([RO.labels(:)], {'trial'}));
trial.RO = RO.data(:,temp);

l = numel(RO.labels); 
for coli = 2:l
    outputcell{1,coli-1} = RO.labels{1,coli};
end


outputcell{1,coli+2} = 'spiketime';
outputcell{1,coli+3} = 'bradi';
outputcell{1,coli+4} = 'normalize';
outputcell{1,coli+5} = '';


% GET ACTIVITY DATA AND NORMALIZE
Adata = RO.data(:,5:end);
d = size(Adata);
Adata = reshape(Adata, [d(1)*d(2) 1]);
AdataN = normalize(Adata, 'range');
AdataN = reshape(AdataN, [d(1) d(2)]);

% MATCH 'spike' TRIALS WITH 'RO' AND GET REARRANGED INDEXES
for ii = 1:size(trial.RO)
rowi = ii+1;
    newidx = find(trial.spike ==  trial.RO(ii));
    
    % LOOP THROUGH THE IMPORTED DATA TO COPY INTO OUTPUT CELL
    for jj = 2:l %loops through the original columns of RO.data
    outputcell{rowi,jj-1} = RO.data(ii,jj);
    end
    
    % ADD NORMALIZED ACTIVITY (at the end) --- non flexible (only for 2
    % columns of data)
    outputcell{rowi, l+4} = AdataN(ii,1); % noramalize activity of dm4
    outputcell{rowi, l+5} = AdataN(ii,2); % noramalize activity of dm2

    % MATCH SPIKE TRIALS 
    if isempty(newidx)
        outputcell{rowi,l+2} = NaN;
        outputcell{rowi,l+3} = NaN;
        
    else
        outputcell{rowi,l+2} = time(newidx);
        outputcell{rowi,l+3} = bradi(newidx);
    end
      
    clear newidx
end

% OUTPUT BRADY DATA REORDERED
source{1, 1} =  'bradi source';
source{1, 2} =  file.spike;

source{2,1} =  'brainA source';
source{2,2} = file.ROmatlab;

source{3,1} = 'normalize: range';

source{5,1} = date;
source{6,1} = 'source'; 
source{6,2} =  mfilename('fullpath');


%% STORE REARRANGED INTO file.ROmatlab
fishref = file.ROmatlab(24:29); %get fish name form standarized name
excelname = fullfile(folder, ['corr03_MATCH_correlation_NormRange' fishref  sheet.ROmatlab '.xls']);

    if isfile(excelname)
        % if file already exists:
        M = strcat ('ALREADY EXISTING FILE!!! / EL ARCHIVO YA EXISTE!!!.', newline , 'The file: "', excelname,  '" already exists in that location.' , newline, 'This code cannot override an existing file.',sprintf("\n"),  'TO SOLVE THE ERROR, change its name, its location or delete before running the code');
        error(M)
    end
    
writecell(outputcell, excelname)
writecell(source, excelname, 'sheet', 'source')


%% Update record 
% 21/11/22: add normalization of the activity data
% Created:
