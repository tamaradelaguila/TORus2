%% PREPARATORY STEPS

% ROOT FOLDER content: 
% pipeline files, 'savecode'
% 'functions' folder
% 'data_structures' folder - for VSDI structures, struct_list
% 'data_bigstructures' folder -for movies and big mat files
% 'plots' folder

%load prespecify settings (change path.rootpath for each computer)
clear
user_settings

%% 1. CREATE LIST OF FISH
% It will be needed for looping across fish
% It has to match the name of the .mat that cointain the 

% load(fullfile(path.data,'grouplist.mat'),'grouplist');
 
grouplist{1} = 'TORus_210320' ;
grouplist{2} = 'TORus_210323' ;
grouplist{3} = 'TORus_210325' ;
grouplist{4} = 'TORus_210409' ;
grouplist{5} = 'TORus_210411' ;
grouplist{6} = 'TORus_210412' ;
grouplist{7} = 'TORus_210421' ;
grouplist{8} = 'TORus_210430' ;
grouplist{9} = 'TORus_210508' ;
grouplist{10} = 'TORus_210509' ;
grouplist{11} = 'TORus_210521' ;
grouplist{12} = 'TORus_210522' ;
save(fullfile(path.data,'grouplist.mat'),'grouplist');

% load(fullfile(path.data,'grouplist.mat'),'grouplist');
%% 2. CREATE VSDI structure (FOR EACH FISH) 
VSDI.ref = 210522 ; %@ SET
VSDI.info.stime = 6; %ms (sampling time) @ SET
VSDI.info.Sonset = 600; % ms (from start of recording) %@ SET
% VSDI.info.Sonset = 590; % ms (from start of recording) %@ SET

% IMPORT LIST. Have to be saved from Brainvision. See notes for details 
listpath =  path.list; %@ SET in user_settings
listpath = fullfile(listpath,strcat('filelist',num2str(VSDI.ref),'.csv'));

list_table=  readtable(listpath);
VSDI.list = table2struct(list_table);
for triali = 1:length(VSDI.list)
VSDI.trialref(triali,1) = str2num(VSDI.list(triali).Name(1:end-5)); %save references for the trials
end

%Get the preceeding ITI for each trial
[VSDI.iti] = get_iti(VSDI); 
TORus('save',VSDI);

%% @ SET
% ADD MANUALLY THE CONDITION FOR EACH TRIAL (in new fields -name them to be able to copy them)...
% Non-included trials: NaN

% Set all to NaN and later add the conditions
for ii = 1:length(VSDI.list)
    VSDI.list(ii).code= NaN;
    VSDI.list(ii).stim_type= NaN;
    VSDI.list(ii).intraderm= NaN;
    VSDI.list(ii).mA= NaN;
    VSDI.list(ii).Sdur= NaN;
  
end

        for ii =1: length(VSDI.list)


              
            if strcmp(VSDI.list(ii).Comment, '0')
              VSDI.list(ii).code = 400;
              VSDI.list(ii).stim_type = 4 ;
              VSDI.list(ii).intraderm = 0 ;

              VSDI.list(ii).mA = 0;

              VSDI.list(ii).Sdur = 0 ;
              
            elseif strcmp(VSDI.list(ii).Comment, '1')
              VSDI.list(ii).code = 401;
              VSDI.list(ii).stim_type = 5 ;
              VSDI.list(ii).intraderm = 0 ;

              VSDI.list(ii).mA = 0.38;

              VSDI.list(ii).Sdur = 0.15 ;

            elseif strcmp(VSDI.list(ii).Comment, '2')
              VSDI.list(ii).code = 402;
              VSDI.list(ii).stim_type = 4 ;
              VSDI.list(ii).intraderm = 0 ;

              VSDI.list(ii).mA = 0.4;

              VSDI.list(ii).Sdur = 0.15 ;
            
            elseif strcmp(VSDI.list(ii).Comment, '3')
              VSDI.list(ii).code = 403;
              VSDI.list(ii).stim_type = 4 ;
              VSDI.list(ii).intraderm = 0 ;

              VSDI.list(ii).mA = 0.6;

              VSDI.list(ii).Sdur = 0.15 ;
              
            elseif strcmp(VSDI.list(ii).Comment, '4')
              VSDI.list(ii).code = 404;
              VSDI.list(ii).stim_type = 4 ;
              VSDI.list(ii).intraderm = 0 ;

              VSDI.list(ii).mA = 0.8;

              VSDI.list(ii).Sdur = 0.15 ;

           end

        end
% 

%AND THEN COPY IT INTO NEW FIELDS (as many as condition - columns you have)

for triali = 1:length(VSDI.list)
VSDI.condition(triali,1) = VSDI.list(triali).code; %@ SET the name of the field so it c an be c opied
VSDI.condition(triali,2) = VSDI.list(triali).stim_type; %@ SET the name of the field so it can be c opied
VSDI.condition(triali,3) = VSDI.list(triali).intraderm; %@ SET the name of the field so it can be c opied
VSDI.condition(triali,4) = VSDI.list(triali).mA; %@ SET the name of the field so it can be c opied
end

VSDI.nonanidx= find(~isnan(VSDI.condition(:,1))) ;

% Get BV times from 'Date' info
for ii = 1:length(VSDI.list)
    date = VSDI.list(ii).Date; %get from VSDI structure
    date.Format = 'HH:mm:ss'; % capital 'H' so it's in 24h format
    date= cellstr(date); %turn into string to keep only hours
    hour = duration(date, 'Format','hh:mm:ss'); % turn again into duration
   
    %store in variable
    VSDI.trialtime(ii).trialref = VSDI.trialref(ii);
    VSDI.trialtime(ii).hour = date;
    VSDI.trialtime(ii).trialref = VSDI.trialref(ii);
end

TORus('save',VSDI)


% VSDI.info.Sdur = 200; % ms(duration of stimulus) %@ SET
for triali = 1:length(VSDI.list)
VSDI.info.Sdur(triali,1) = VSDI.list(triali).Sdur; %@ SET the name of the field so it c an be c opied
end

% Save changes
TORus('save',VSDI);

% test saving
clear
[VSDI] = TORus('load',8);
