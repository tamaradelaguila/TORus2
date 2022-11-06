%% Load pike matfile
% Previous steps in 'spike_export.txt'

nfish = 5; 
user_settings
VSDI = ROCus('load',nfish);


% ----------------------------------------------------
%% get reference file with [filename initial-time] in each row
% ----------------------------------------------------
spike.ref = VSDI.ref;

pathspike= '/home/tamara/Documents/MATLAB/VSDI/ROCtone/data/dataspike'; 

% spikeref = readtable(fullfile(pathspike,num2str(VSDI.ref),[num2str(VSDI.ref) 'spikelist']));
file =fullfile(pathspike,'210920','210920spikelist.csv');

spikeref=  readtable(file);
spikeref=  table2cell(spikeref);


% GET REFERENCE TIME FOR LATER 
spike.reftime.trial = spikeref{1,1};%trial to reference from
spike.reftime.spiketime = spikeref{1,2}; %time frmo the from spike (setting cursor)

% ----------------------------------------------------
% IMPORT VALUES FROM ALL SPIKE FILES AND CONCATENATE

spike.ro = [];

spike.ecg = [];
spike.ecg_timebase = [];

for ii = 2:length(spikeref)
    % load spike file
%     filepath = fullfile(pathspike, num2str(VSDI.ref), spikeref{ii,1});
        filepath = fullfile(pathspike, '210920', spikeref{ii,1});
        t0 = spikeref{ii,2};
        load(filepath)
        
        spike.ro = cat (1, spike.ro, t0 + seconds(ro.times));
        
%         spike.ecg= cat(1, spike.ecg, downsample (ecg.values, 10));
%         temp_times = downsample(ecg.times,10) ;
%         spike.ecg_timebase = cat(spike.ecg_timebase, t0+temp_times);
%         
%         spike.cs5 = ;
%         spike.cs1 = ;
%         spike.us = ;
%         
clear temp_times
end


%% MATCH - GET REFERENCE TIMES AND ESTIMATED ROtimes correspondence

% MANUALLY ADD REFERENCE

spike.reftime.trialidx =  find(VSDI.trialref == sscanf(spike.reftime.trial, '%d_%d') ); %find the corresponding matlab index
bvtime=  VSDI.list(spike.reftime.trialidx).Date; %from VSDI.list (BVfile time of creation)
bvtime.Format = 'HH:mm:ss';
spike.reftime.BVtime=bvtime;

% gap between 
gap = spike.reftime.BVtime - spike.reftime.spiketime;

for ii= 1:length(ro.times)
spike.match(ii).ROspike = seconds(ro.times(ii)); % absolute seconds

spike.match(ii).BVestim = spike.match(ii).ROspike + spike.reftime.gap;

spike.match(ii).ROspike_h = spike.reftime.spike0 + spike.match(ii).ROspike; %
spike.match(ii).BVestim_h = spike.reftime.spike0 + spike.match(ii).BVestim; 

end

spike.reftime.gap =gap; 

%% MATCH - GET timeinfo from BV files

% 1. GET Original BV files time of creation
%     for ii= 1:length(VSDI.trialtime)
%     BVtimes_sec(ii) = duration(VSDI.trialtime(ii).hour, 'Format', 'hh:mm:ss')- spike.reftime.spike0; 
%     end
%     BVtimes_sec = BVtimes_sec';
%     BVtimes_sec.Format = 's';
    
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
    
%% MATCH - MATCHING PROCESS

% ------------------------------------------------
% APPLY MATCHING FUNCTION AND SAVE CORRESPONDENCE
% ------------------------------------------------

% For each spike-ROtime, find the closest 
for ii= 1:length(spike.match)
idx= find_closest_timeidx(spike.match(ii).BVestim, BVtimes_sec);

spike.match(ii).BVtrialidx =  idx;
spike.match(ii).BVfiletime = BVtimes_sec(idx);
spike.match(ii).BVdml = VSDI.trialref(idx);
end

% Turn into string to easily read 
output{1,1}= 'ROspike_h';
output{1,2}= 'BVestim_h';
output{1,3}= 'trialidx';
output{1,4}= 'dml file';
output{1,5}= 'manual correction';

for ii=1:length(spike.match)
output{ii+1,1}= spike.match(ii).ROspike_h;
output{ii+1,2}= spike.match(ii).BVestim_h;
idx = spike.match(ii).BVtrialidx;
output{ii+1,3}=  idx; 
output{ii+1,4}=  VSDI.trialref(idx);
end

% check whether there is a repeated matched trial:
            for ti = 2:length(spike.match)
                if spike.match(ti).BVtrialidx- spike.match(ti-1).BVtrialidx < 1
                    output{ti-1,6}= 'check'; 
                end
            end

%%  FEEDBACK
excelfile =fullfile(pathspike,'210920','test.xls');
writecell(output, file, 'Sheet', 'matlab')

% ------------------------------------------------------------------------
% MANUALLY CHECK AND CORRECT IN THE OUTPUT EXCEL (correct in the column
% 'manual_correction')
% ------------------------------------------------------------------------

checked_excel=  xlsread(excelfile, 'checked');

spike.match_ok.ROspike_h = spike.match.ROspike_h;
spike.match_ok.BVdml = spike.match.BVdml;
spike.match_ok.trialidx = spike.match.BVtrialidx;

% substitute the corrected values
excelidx = find(~isnan(checked_excel(:,3)));
for ii = makeRow(excelidx)
    oldmatch = (checked_excel(ii,2);
    newmatch = checked_excel(ii,3);
    % find the index in the match list to replace
    listidx = find(spike.match_ok.BVdml == oldmatch);
    
    spike.match_ok(listidx).BVdml = newmatch ; %substitute the match in the list

    % find the matlab index corresponding to the new matched dml file 
    newidx = find(VSDI.trialref == newmatch); 
    spike.match_ok(listidx).trialidx = newidx;
    
end


%% TRANSFER INTO VSDI.structure so the spike times can be used

for ii = 1:length(spike.match)
trialidx = spike.match_ok(ii).trialidx; 
VSDI.spike.RO(trialidx,1) = round(seconds(spike.match_ok(ii).ROspike_h),2); %save in the correspondent BVindex
VSDI.spike.BVdml(trialidx,1) = spike.match_ok(ii).BVdml; % 
end
% VSDI matching structure


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

