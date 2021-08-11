%% ABSOLUTE THRESHOLD FROM ROI WAVE
% If any roi wave (from selected ones) is above the threshold

abs_thr = 0.4; 
reject_idx = [];

for triali = VSDI.nonanidx
    if max(wave(:,triali)) > abs_thr 
    reject_idx = [reject_idx; triali];
    end
end

fieldname = ['abs',num2str(abs_thr)];
VSDI.reject.(fieldname)= reject_idx;

%% ABSOLUTE THRESHOLD ON GS rejection criterion
% If the GS is above the threshold
clear 
user_settings 
pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/'; 

% nfish = TORus('nsubject', 210409)
for nfish = [2 3 4 8]
    VSDI = TORus('load', nfish);
VSDroiTS = TORus('loadwave',nfish);
abs_thr = 0.25; 
reject_idx = [];

wave = VSDroiTS.filt306.GS; 

for triali = makeRow(VSDI.nonanidx)
    if max(wave(:,triali)) > abs_thr
    reject_idx = [reject_idx; triali];
    end
end

fieldname = ['GSabs',num2str(abs_thr)]; fieldname = strrep(fieldname,'.','');
VSDI.reject.(fieldname)= reject_idx;


% plot rejected: 
plot_rejected(wave, makeRow(VSDI.nonanidx), VSDI.reject.(fieldname))
title([num2str(VSDI.ref),'rejected by criterion:', fieldname])

saveas(gcf, fullfile(pathplot, [num2str(VSDI.ref), 'rejected', fieldname, '.jpg']),'jpg')

close 
TORus('save', VSDI);

end
%% CUMMULATIVE ERROR DEVIATION FRMO MEAN IN GS rejection criterion
% GLOBAL SIGNAL METHOD
%  BASED ON THE AVERAGE SIGNAL OF THE WHOLE CROPPED BRAIN IN THE TRIALS
%  INCLUDED IN THE STUDY (VSDI.trials_in)
% Method by (Cheng, 2006) applied to the whole brain.

% Extract whole brain averaged (GS) from all trials from each conditions
% and average across trials
% substract from 

% STEP 0 Preprocess and crop brain (p10)
% STEP 1 GS for each trial. Average across trials (totalGS)

% *function to compute totalGS and total-Std.
% STEP 2 (for each trial): extract totalGS, compute the standart deviation of the residuals
% STEP 3 : compute the mean of the SD from each trial
% STEP 4 (for each trial): Discard if std of the residuals for that trial
% is > 2 times the mean of the SD
%^* in a loop 

% If the GS is above the threshold

std_factor = 2;
saveplot = 1; 
pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/'; 

for nfish = [2 3 4 8]
VSDI = TORus('load',nfish);
VSDroiTS = TORus('loadwave',nfish);
onsetidx = find(VSDI.timebase == 0);

%copy GS
GS = VSDroiTS.filt306.GS;

% ---------------------------------------------------------
% FOR BASELINE --------------------------------------------
% find aberrant trials in baseline period for all trials


all_idx = VSDI.nonanidx; % it will include all trials from all conditions

[out_basel] = find_aberrant(GS(1:onsetidx,:), all_idx,  std_factor, '0'); title([num2str(VSDI.ref) '.Basel GS from rejected(grey)from all included trials'])

        if saveplot == 1
           
                name = strcat(num2str(VSDI.ref),'GSbasel_rejected'); 
                saveas(gcf, fullfile(pathplot, [name, '.jpg']),'jpg')
                close all
        end

% -------------------------------------------------------
% FOR POST-S --------------------------------------------

% find aberrant trials for each condition in the post-Stimulus period
cond_codes = unique(VSDI.condition(:,1)); cond_codes = cond_codes(~isnan(cond_codes)); 

out_specific = [];

for ci = 1:numel(cond_codes)
codi  = cond_codes(ci);
idx = find(VSDI.condition(:,1) == codi); 
temp_aberrant = find_aberrant(GS, idx, std_factor, 1); title([num2str(VSDI.ref) '.Post-S GS  from rejected(grey)from cond' num2str(cond_codes(ci))])
out_cond = [out_specific temp_aberrant];
end

        if saveplot == 1
           
                h =  findobj('type','figure');
                n = length(h);

                for ii = 1:n
                name = strcat(num2str(VSDI.ref),'GSpost_rejected', 'c', num2str(h(ii).Number)); 
                saveas(h(ii), fullfile(pathplot, [name, '.jpg']),'jpg')
                close(h(ii))
                end
                clear h n
 
        end


% store in structure 
fieldn = ['GSdeviat', num2str(std_factor),'sd']; fieldn= strrep(fieldn,'.','');
VSDI.reject.(fieldn) = sort(union(out_basel, out_specific));

TORus('save', VSDI)

clear   aberrant_idx reject rejectB rejectP...
        aberrant_idxall aberrant_idx...
        GS VSDI 
end
%% PLOT ALL REJECTED from certain conditions

pathplot = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/'; 

for nfish =  [2 3 4 8]
VSDI = TORus('load', nfish);
VSDroiTS = TORus('loadwave', nfish);
waves = VSDroiTS.filt306.GS;
totalout = sort(unique([VSDI.reject.GSabs025 VSDI.reject.GSdeviat2sd]));

plot_rejected(waves) ; title([num2str(VSDI.ref), '- GS of rejected (GSdeviat2.5std + abs0.25 thresh)'])
    saveas(gcf, fullfile(pathplot, [num2str(VSDI.ref), 'rejected', fieldname, '.jpg']),'jpg')
    close 
end 

%% PRINT REJECTED IN EXCEL

clear

for nfish =  [2 3 4 8]
VSDI = TORus('load',nfish);

out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject';
out.name = 'rejected';
export.excel = 1;

rejected = sort(unique([VSDI.reject.GSabs025 ; VSDI.reject.GSdeviat2sd])) ;

localoutput(:,3) = rejected;
localoutput(1:end,2) =NaN;

localoutput(:,1) = VSDI.trialref(rejected);
localoutput(:,4) = NaN;
localoutput(:,5) = VSDI.condition(rejected,1); %condition

    % WRITE EXCEL
    if export.excel == 1
        out.sheet = [num2str(VSDI.ref)];
    excelname = fullfile(out.folder,['rejected_trials',out.name,'.xlsx']);
    writematrix (localoutput,excelname,'sheet',out.sheet)
    end
    clear localoutput rejected
end %nfish

%% PLOT REJECTED FOR EACH CONDITION AND SAVE
    
clear
user_settings

for nfish =  [2 3 4 8]
[VSDI] = TORus('load',nfish);

temp = TORus('loadwave',nfish');
wavesGS  = temp.filt306.GS;

% windows in which analyse the measures

cond_codes = unique(VSDI.condition(:,1));
cond_codes = cond_codes(~isnan(cond_codes));

outidx = sort(unique([VSDI.reject.GSabs025; VSDI.reject.GSdeviat2sd]));

for ci = 1:numel(cond_codes)
    sel_trials = find(VSDI.condition(:,1) == cond_codes(ci));
    plot_rejected(wavesGS, sel_trials, outidx, [-0.1 0.4])
    title([num2str(VSDI.ref),': total of rejected (red) cond.', num2str(cond_codes(ci)) ])
%     
    out.name = [num2str(VSDI.ref) '_plot_reject_cond',num2str(cond_codes(ci))];
    pathsaveR = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/reject/'; 
    saveas(gcf, fullfile(pathsaveR, out.name), 'jpg')
    close 
    clear seltrials 
end 
clear outidx temp waveGS
end
blob()