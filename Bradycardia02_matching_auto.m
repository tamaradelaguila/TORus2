%% Load spike matfile
% Previous steps in 'spike_export.txt'
% Bradycardia can be previously calculated from the spike file with the
% code BRADYCARDIA_flexible_counting_code.m


clear
nfish = 13; 
user_settings

% Spike matfile to process stims
spikefile = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike/220611spike.mat' ;

% Csv with brady data
sourcebrady.in = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike/excels_brady'; 
sourcebrady.filelist= '220611spike_BRADYcodeinput.csv';

VSDI = TORus('load',nfish);

%% MATCH STEP 1- GET REFERENCE TIMES AND ESTIMATED ROtimes correspondence
% ....................................................................

% MANUALLY ADD REFERENCE
spike.reftime.trial = '007A';%trial to reference from
spike.reftime.spiketime = '14:54:20'; % corresponing spike time (setting cursor on the RO onset) - ATT: time of the day
spike.reftime.spike0 = '14:51:19';  %starting time from spike 
spike.reftime.trialidx = 8 ;%corresponding BrainVision index for that trial

% % (alternative to manually add: automatical finding of the corresponding
% % VSDI.trialref index from spike.reftime.trial
% % spike.reftime.trialidx =  find(VSDI.trialref == sscanf(spike.reftime.trial, '%d_%d') ); %find the corresponding matlab index
% temp = VSDI.ref*1000 + sscanf(spike.reftime.trial, '%d_%d'); %conversion to the format of the trial reference
% spike.reftime.trialidx =  find(VSDI.trialref == temp ); %find the corresponding matlab index

% ... and turn into 'duration' vectors
bvtime=  VSDI.list(spike.reftime.trialidx).Date; %from VSDI.list (BVfile time of creation)
% bvtime.Format = 'HH:mm:ss';
bvtime =  datestr(bvtime, 'HH:MM:SS');
spike.reftime.BVtime=duration(bvtime, 'inputformat', 'hh:mm:ss');
clear bvtime

spike.reftime.spiketime = duration(spike.reftime.spiketime, 'Format','hh:mm:ss');
spike.reftime.spike0 = duration(spike.reftime.spike0, 'Format','hh:mm:ss');

% CALCULATE GAP BETWEEN COMPUTERS 
% gap = spike.reftime.BVtime - spike.reftime.spiketime;
spike.reftime.gap =   spike.reftime.spiketime - spike.reftime.BVtime;

% LOAD THE .csv FILE WITH THE BRADICARDIA 
% Previously extracted from: /home/tamara/Documents/MATLAB/CC_brady_counting_codes/BRADYCARDIA_flexible_counting_code.m
spiketemp = load(spikefile);
ROtimes = spiketemp.RO.times;

% the first row should the the seconds of arrival of the stim
% Assign the variables depending on which measure was used:

% Get the 
for ii= 1:length(ROtimes)
spike.match(ii).ROspike = seconds(ROtimes(ii)); % absolute seconds

spike.match(ii).BVestim = spike.match(ii).ROspike + spike.reftime.gap;

spike.match(ii).ROspike_h = spike.reftime.spike0 + spike.match(ii).ROspike; %
spike.match(ii).BVestim_h = spike.reftime.spike0 + spike.match(ii).BVestim; 

end

% Save as structure to easily visualize:
spike.matchtable = struct2table(spike.match);
%writetable(struct2table(spike.match), 'Structure_Example.csv')

%% MATCH STEP 2 - GET timeinfo from BV files
% ....................................................................

% 1. GET Original BV files time of creation
    for ii= 1:length(VSDI.trialtime)
    BVtimes_sec(ii) = duration(VSDI.trialtime(ii).hour, 'Format', 'hh:mm:ss')- spike.reftime.spike0; 
    end
    BVtimes_sec = BVtimes_sec';
    BVtimes_sec.Format = 's';
    
% 2. MANUALLY CORRECT: BEFORE MATCHING CHECK FOR INCONGRUENCIES and correct them into the VSDI.trialtime:
    ...check if there is any odd 'iti', or any 'BVtime' that violates the flow of time (i.e., that is earlier than the previous one), and correct the 'time' mismatch from the time in the spike file
    ... skip the non-interest trials (that is, non-recorded trials)
    ... add a new colum/field with 'old_hour' and correct it in the first 'hour' one

% check if there is a missing BV file  
if length(find(diff(VSDI.trialref)==0)) > 0
	disp('"there is at least 1 missing BV file"')

else 
    disp('"there is no missing BV file"')
    
end   

% check if there is any 
    check_idx = [];
    for ti = 2:length(BVtimes_sec)
  
        if VSDI.list(ti).Date < VSDI.list(ti-1).Date
            check_idx = [check_idx; ti-1];
            spike.match(ti-1).check= 'check odd iti'; 
        end
    end
disp(['CHECK INDEXES:' num2str(makeRow(check_idx))])

    if ~isempty(check_idx)
        disp(['"check iti in trials idx:"' makeRow(check_idx)])
    else
        disp('"there is no odd iti"')
    end
    

% % 3. AND SUBSTITUTE THE CORRECT VALUES INTO THE BVtimes_sec
% 
%     % GET timeinfo from BV files
%     for ii= 1:length(VSDI.trialtime)
%     BVtimes_sec(ii) = duration(VSDI.trialtime(ii).hour, 'Format', 'hh:mm:ss')- spike.reftime.spike0;
%     end
%     BVtimes_sec = BVtimes_sec';
%     BVtimes_sec.Format = 's';
    
    ...and run again to see if some mismatch remain
            check_idx = [];
            for ti = 2:length(BVtimes_sec)
                if seconds(BVtimes_sec(ti)) < seconds(BVtimes_sec(ti-1))
                    check_idx = [check_idx; ti-1];
                end
            end
disp(['...Still need to check indexes:' num2str(makeRow(check_idx))])


%% MATCH STEP 3 (OPTION 1 - FROM CSV ) - MATCHING PROCESS ITSELF
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
VSDI.brady.data = date; % UPDATED DATE WHEN IT WAS ADDED
% TORus('save', VSDI)


%% TO FURTHER USE THE BRADY MEASURES

% To find a specific trial measure: first find the idx of the bradydata
idx = find([VSDI.brady.data(:).idx] == 33);

% loop through bradylist
for ii= 1:length(VSDI.brady.measurelist) 
field =  bradylist{ii};
bradi = VSDI.brady.data.(field); 
end



%% MATCH STEP 3 (OPTION 2 AUTOMATIC) - MATCHING PROCESS ITSELF
% % ....................................................................
% % !!! CHECK PROPERLY: This code only works when the BVtrials are saved short after acquired
% FINISH TO MATCH THE OUTPUT WITH OPTION 1

% % ------------------------------------------------
% % APPLY MATCHING FUNCTION AND SAVE CORRESPONDENCE
% % ------------------------------------------------
% 
% % For each spike-ROtime, find the closest 
% for ii= 1:length(spike.match)
% idx= find_closest_timeidx(spike.match(ii).BVestim, BVtimes_sec);
% 
% spike.match(ii).BVfiletime = BVtimes_sec(idx);
% spike.match(ii).BVdml = VSDI.trialref(idx);
% spike.match(ii).BVtrialidx =  idx;
% 
% end
% spike.matchtable = struct2table(spike.match);
% 
% % Check whether there is a repeated matched trial:
%             for ti = 2:height(spike.matchtable)
%                 if spike.match(ti).BVtrialidx- spike.match(ti-1).BVtrialidx < 1
%                     spike.matchtable{ti-1,8}= {'check'}; 
% 
%                 end
%             end
% % ------------------------------------------------------------------------
% % ------------------------------------------------------------------------
% %  FEEDBACK
% excelfile =fullfile(pathspike,'210920','test.xls');
% writecell(output, file, 'Sheet', 'matlab')
% 
% % ------------------------------------------------------------------------
% % MANUALLY CHECK AND CORRECT IN THE OUTPUT EXCEL (correct in the column
% % 'manual_correction')
% % ------------------------------------------------------------------------
% 
% checked_excel=  xlsread(excelfile, 'checked');
% 
% spike.match_ok.ROspike_h = spike.match.ROspike_h;
% spike.match_ok.BVdml = spike.match.BVdml;
% spike.match_ok.trialidx = spike.match.BVtrialidx;
% 
% % substitute the corrected values
% excelidx = find(~isnan(checked_excel(:,3)));
% for ii = makeRow(excelidx)
%     oldmatch = checked_excel(ii,2);
%     newmatch = checked_excel(ii,3);
%     % find the index in the match list to replace
%     listidx = find(spike.match_ok.BVdml == oldmatch);
%     
%     spike.match_ok(listidx).BVdml = newmatch ; %substitute the match in the list
% 
%     % find the matlab index corresponding to the new matched dml file 
%     newidx = find(VSDI.trialref == newmatch); 
%     spike.match_ok(listidx).trialidx = newidx;
%     
% end

% OUTPUT TO CHECK
% tempname = ['CHECK_MATCH_' num2str(VSDI.ref)];
% writetable(struct2table(spike.match_ok), fullfile(sourcebrady.in,tempname), 'filetype', 'spreadsheet' )


%% CHECKING - TEST WHETHER THE PRESENCE OF STIMULUS IN SPIKE COINCIDES WITH THE PRESENCE OF STIMULUS IN BRAINVISION
% And whether there are spike trials not recorded on BV
% 
% check_trialidx=[];
% for ii= 1:makeRow(VSDI.nonanidx)
%     % check if there is stimulus in that spike trial
%     
%     limit1 = spike.Sonset >= VSDI.spike.RO(ii);
%     coinc = spike.Sonset== VSDI.spike.RO(ii);
%     limit2 = spike.Sonset < VSDI.spike.RO(ii)+2;
%     
%     spike_idx = find (( limit1 & limit2 )); 
%     if ~isempty(spike_idx) && spike_idx > 0
%     stim_spike_flag = 1;
%     else 
%     stim_spike_flag = 0;
%     end
%     
%     % check if there is stimulus according to BrainVision
%          stim_BV_flag = VSDI.list(ii).mA ~=0 && ~isnan(VSDI.list(ii).mA );
%         
%     % if they do not coincide: check index
% if stim_spike_flag ~= stim_BV_flag
%     check_trialidx=[check_trialidx; ii];
% end
% 
% clear stim_spike_flag stim_BV_flag
% end

%% SAVE SPIKE STRUCTURE
% TORus('save',VSDI) %DEMUTE
% TORus('savespike',spike)

%% COUNTING SPIKES 

% PLOT - individual trial tiles + waves + heart example >> do the three
% fish

% first iti(s) after stimulus onset

%% NOTE THAT ONE OF THE FOLLOWING ISSUES MIGHT HAPPEN:

% that the BVfile is saved just before or after the following trial
% that there can be a spike event without a real BV file (but with no
% missing BVfile, because it was rewritten)
% that there can be any missing BV file

%% DEPRECATED:
% 
% % Deprecated because there will be one brady file per spike file.
% % If concatenation is needed, it will have to be handled differently (for
% % instance, concatenate the bradylists before processing
% % --------------------------------------------------------------
% % CONCATENATE IF NEEDED (if there are multiple files)
% % --------------------------------------------------------------
% % get reference file with [filename initial-time] in each row
% 
% spike.ref = VSDI.ref;
% 
% pathspike= '/home/tamara/Documents/MATLAB/VSDI/ROCtone/data/dataspike'; 
% 
% % spikeref = readtable(fullfile(pathspike,num2str(VSDI.ref),[num2str(VSDI.ref) 'spikelist']));
% file =fullfile(pathspike,'220611','220611spikelist..csv');
% 
% spikeref=  readtable(file);
% spikeref=  table2cell(spikeref);
% 
% % GET REFERENCE TIME FOR LATER 
% spike.reftime.trial = spikeref{1,1};%trial to reference from
% spike.reftime.spiketime = spikeref{1,2}; %time from spike (setting cursor)
% 
% % ----------------------------------------------------
% % IMPORT VALUES FROM ALL SPIKE FILES AND CONCATENATE
% 
% spike.ro = [];
% 
% spike.ecg = [];
% spike.ecg_timebase = [];
% 
% for ii = 2:length(spikeref)
%     % load spike file
% %     filepath = fullfile(pathspike, num2str(VSDI.ref), spikeref{ii,1});
%         filepath = fullfile(pathspike, '210920', spikeref{ii,1});
%         t0 = spikeref{ii,2};
%         load(filepath)
%         
%         spike.ro = cat (1, spike.ro, t0 + seconds(ro.times));
%         
% %         spike.ecg= cat(1, spike.ecg, downsample (ecg.values, 10));
% %         temp_times = downsample(ecg.times,10) ;
% %         spike.ecg_timebase = cat(spike.ecg_timebase, t0+temp_times);
% %         
% %         spike.cs5 = ;
% %         spike.cs1 = ;
% %         spike.us = ;
% %         
% clear temp_times
% end

%% UPDATE HISTORY

% 24/08/22 Finish code for Matching option 2 - manual 
% 19/08/22  Created from: spike_processing_easyfeedback
