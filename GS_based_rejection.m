%% TRIAL REJECTION 

%% GLOBAL SIGNAL METHOD
%  BASED ON THE AVERAGE SIGNAL OF THE WHOLE CROPPED BRAIN IN THE TRIALS
%  INCLUDED IN THE STUDY (VSDI.trials_in)
% Method by (Cheng, 2006) applied to the whole brain.

% Extract whole brain averaged (GS) from all trials from each conditions
% and average across trials
% substract from 

% STEP 0 Preprocess and make cropmask (s01), extract GS (s04)
% STEP 1 GS for each trial. Average across trials (totalGS)

% *function to compute totalGS and total-Std.
% STEP 2 (for each trial): extract totalGS, compute the standart deviation of the residuals
% STEP 3 : compute the mean of the SD from each trial
% STEP 4 (for each trial): Discard if std of the residuals for that trial
% is > 2 times the mean of the SD
%^* in a loop 

clear
user_settings

saveplot = 1; % of rejected trials

for nfish =1
VSDI = TORus('load',nfish);
VSDroiTS= TORus('loadwave',nfish);
% load VSDI and VSDroiTS structure (to extract GS)

% ---------------------------------------------------------
% FOR BASELINE --------------------------------------------
% find aberrant trials in baseline period for all trials
all_idx = VSDI.nonanidx;
idx_onset = find();

[reject_baseline] = find_aberrant(GS(1:,:), all_idx,  1); title([struct_list{nfish}(5:end-4) '.Basel GS from rejected(grey)from all trials'])

        if saveplot == 1
           
                name = strcat('GSbasel_rejected', num2str(VSDI.ref)); 
                saveas(gcf, fullfile(pathplot, [name, '.jpg']),'jpg')
                close all
        end

% -------------------------------------------------------
% FOR POST-S --------------------------------------------

% find aberrant trials for each condition in the post-Stimulus period
c0idx = intersect(find(VSDI.condition==0), VSDI.trials_in);
c1idx =  intersect(find(VSDI.condition==1), VSDI.trials_in);
c2idx =  intersect(find(VSDI.condition==2), VSDI.trials_in);

[aberrant_idx0p] = find_aberrant(GS(201:end-2,:), c0idx, 1); title([struct_list{nfish}(5:end-4) '.Post-S GS  from rejected(grey)from cond 0'])
[aberrant_idx1p] = find_aberrant(GS(201:end-2,:), c1idx, 1); title([struct_list{nfish}(5:end-4) '.Post-S GS postSfrom rejected(grey)from cond 1'])
[aberrant_idx2p] = find_aberrant(GS(201:end-2,:), c2idx, 1); title([struct_list{nfish}(5:end-4) '.Post-S GS from rejected(grey)from cond 2'])

        if saveplot == 1
           
                h =  findobj('type','figure');
                n = length(h);

                for ii = 1:n
                name = strcat('GSpost_rejected', num2str(VSDI.ref),'c', num2str(h(ii).Number)); 
                saveas(h(ii), fullfile(pathplot, [name, '.jpg']),'jpg')
                close(h(ii))
                end
                clear h n
 
        end

rejectP = sort([aberrant_idx0p aberrant_idx1p aberrant_idx2p]);

% store in structure 
VSDI.reject.GStotal = sort(union(rejectB, rejectP)); 
VSDI.reject.GSbasel= rejectB;
VSDI.reject.GSpost = rejectP;

save (strcat(VSDI.pathsave,num2str(VSDI.ref),'.mat'), 'VSDI')

clear c0idx c1idx c2idx aberrant_idx0 aberrant_idx1 aberrant_idx2 reject rejectB rejectP...
    aberrant_idxall aberrant_idx0p aberrant_idx1p aberrant_idx2p
clear GS VSDI 
end

%% PEAKS FROM GS METHOD (SHARK-LIKE rejection)

for nfish =2:7
    rootpath = 'C:\Users\User\Documents\UGent_brugge\';
    pathdrive= 'A:\TOR_erps\ANALYSIS\data_bigstructures';

    load (fullfile(rootpath, 'data_structures','struct_list.mat' ))
VSDI = loadfish(nfish, rootpath);

% load VSDI structure and movie (to extract GS)
% load(fullfile(pathbig, ['moviep10_' struct_list{nfish}(5:end)]))
load(fullfile(pathdrive, ['moviep10_' struct_list{nfish}(5:end)]))
%copy GS
GS = VSDmoviep10.GS(1:680,:); clear VSDmoviep10; 

saveplot =1; 
pathplot = fullfile(rootpath, 'plots', 'rejected_sharks');
% FIND SHARKS (it does it does it for baseline and condition-specific
 % directly)  
% SDthresh = 2;
% sharks =  find_shark (GS, VSDI, SDthresh, 1);
%         if saveplot == 1
%                 title([num2str(VSDI.ref), 'sharks. SD=', num2str(SDthresh)])
%                 name = strcat('sharks_rejected', num2str(VSDI.ref)); 
%                 saveas(gcf, fullfile(pathplot, [name, '.jpg']),'jpg')
%                 close all
%         end

%  FIND TRIALS THAT CROSS AN ABSOLUTE THRESHOLD    
abs_thresh =  0.35;
abovethreshold = find_suprathresh(GS, VSDI, abs_thresh, 1);

        if saveplot == 1
                title([num2str(VSDI.ref), 'Above threshold=', num2str(abs_thresh)])
                name = strcat('Reject_abovethresh', num2str(VSDI.ref), '_thresh',num2str(abs_thresh)); 
                saveas(gcf, fullfile(pathplot, [name, '.jpg']),'jpg')
                close all
        end


% store in structure 
% VSDI.reject.sharks = sharks;
VSDI.reject.absthresh = abovethreshold;
save (strcat(VSDI.pathsave,num2str(VSDI.ref),'.mat'), 'VSDI')

clear 
end


for nfish =2:7
    rootpath = 'C:\Users\User\Documents\UGent_brugge\';
    pathdrive= 'A:\TOR_erps\ANALYSIS\data_bigstructure s';

    load (fullfile(rootpath, 'data_structures','struct_list.mat' ))
VSDI = loadfish(nfish, rootpath);

% load VSDI structure and movie (to extract GS)
% load(fullfile(pathbig, ['moviep10_' struct_list{nfish}(5:end)]))
load(fullfile(pathdrive, ['moviep10_' struct_list{nfish}(5:end)]))
%copy GS
GS = VSDmoviep10.GS(1:680,:); clear VSDmoviep10; 



    exclude = union (VSDI.reject.GStotal, VSDI.reject.sharks); exclude = union (exclude, VSDI.reject.absthresh);

    figure
    plot(GS(:,exclude), 'color','#97978F'); hold on;
    plot(GS(:,setdiff(VSDI.trials_in,exclude)), 'color','#FBC4BF'); 

                title([num2str(VSDI.ref), 'All rejected'])
                name = strcat('All_rejected', num2str(VSDI.ref)); 
                saveas(gcf, fullfile(rootpath, 'plots','rejected_sharks', [name, '.jpg']),'jpg')
                
                close 


clear 
end

%% MANUAL REJECTION

%fish 16
VSDI = loadfish(2); 
VSDI.reject.manual = [18 28 31 38 50 51 55 73 76 90 93 94 97 111 112 118];
save(strcat(VSDI.pathsave, num2str(VSDI.ref)),'VSDI')
clear

%fish 17
VSDI = loadfish(3); 
VSDI.reject.manual = [17 19 26 30 40 46 52 86];
save(strcat(VSDI.pathsave, num2str(VSDI.ref)),'VSDI')
clear

%fish 18
VSDI = loadfish(4); 
VSDI.reject.manual = [12 13 16 18 20 23 61 63 65 66 80 96];
save(strcat(VSDI.pathsave, num2str(VSDI.ref)),'VSDI')
clear

%fish 19
VSDI = loadfish(5); 
VSDI.reject.manual = [18 19 63 91 99 101 128];
save(strcat(VSDI.pathsave, num2str(VSDI.ref)),'VSDI')
clear

%fish 20
VSDI = loadfish(6); 
VSDI.reject.manual = [17 35 36 46 50 65 72 92 99 112 120 123 124 128 131 132];
save(strcat(VSDI.pathsave, num2str(VSDI.ref)),'VSDI')
clear

%fish 21
VSDI = loadfish(7); 
VSDI.reject.manual = [14 18 127];
save(strcat(VSDI.pathsave, num2str(VSDI.ref)),'VSDI')
clear

%% SHARK CHARACTERIZATION
clear

% Paths
rootpath = 'C:\Users\User\Documents\UGent_brugge\';
load(fullfile(rootpath, 'data_structures','struct_list.mat'))

pathreg = fullfile(rootpath, 'data_structures');
pathbig = fullfile(rootpath, 'data_bigstructures');

pathdrive= fullfile('A:','TOR_erps','ANALYSIS','data_bigstructures');


% LOAD (1) VSDI, (2) TS and (3) connectivity-matrix structures, (4) movie
% (to extract GS) and (5) imreg information 
nfish = 2;
[VSDI, TSroi] = loadfish(nfish, rootpath);
% load connectivity matrices
namemat= strcat('Fconn1w_', struct_list{nfish});
load(fullfile(rootpath, 'data_structures', namemat)); connmatrix = Fconn_1w.datap10;
%load movie, copy GS and clear big file
load(fullfile(pathdrive, ['moviep10_' struct_list{nfish}(5:end)]))
GS = VSDmoviep10.GS(1:680,:); clear VSDmoviep10;

load(fullfile('data_structures','others','imregform_200616.mat')) %total amount of registration correction respect to the previous frame needed

for ti =  VSDI.trials_in'
sgtitle (['explore sharks. ', num2str(VSDI.ref),'.trial idx', num2str(ti)])

%% 
subplot(2,2,1)
     plot (VSDI.timebase, TSroi.datap10(1:680,:,ti)); hold on
     plot(VSDI.timebase,GS(1:680,ti), 'color','k', 'linewidth', 0.8); hold off %range for 1 window (post-S)
     set(gca, 'xlim', [-1200 2880], 'ylim',[-0.2 0.4])
    xline(0, '--k');
    title('ROI timeser & GS(black)')
    
subplot(2,2,2)
    imagesc(connmatrix(:,:,ti)); colorbar; set(gca, 'clim', [-1 1]); colormap(polarmap(parula))
    title('funct connect (postS)')
    
subplot(2,2,3)
     plot (VSDI.timebase, transf(ti,:));
     set(gca, 'xlim', [-1200 2880])

    title('sum regist respect previous frame')

    subplot(2,2,4)
    meanwind = movsum(transf(ti,:), 13); 

     plot (VSDI.timebase, meanwind);
     set(gca, 'xlim', [-1200 2880], 'ylim', [0 0.25])

    title('sum regist mean ')

% % Size settings in case of 3 connection matrices
%  set(gcf, 'PaperUnits', 'centimeters');
%  x_width=16;y_width=16;
%  set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% 
saveas(gcf,fullfile(rootpath, 'plots','rejection_explore', ['sharks_charact_',num2str(VSDI.ref),'trial_',num2str(ti) ,'.jpg']),'jpg')
close
% pause
end

%% ######################################################################################
% FOR WHOLE TRIAL - DISCARDED CODE for GS method --------------------------------------------
% find aberrant trials for each condition
% [aberrant_idx0] = find_aberrant(GS, c0idx, 1); title([struct_list{nfish}(5:end-4) '.GS from rejected(grey)from cond 0'])
% [aberrant_idx1] = find_aberrant(GS, c1idx, 1); title([struct_list{nfish}(5:end-4) '.GS from rejected(grey)from cond 1'])
% [aberrant_idx2] = find_aberrant(GS, c2idx, 1); title([struct_list{nfish}(5:end-4) '.GS from rejected(grey)from cond 2'])
%         if saveplot == 1
%                 h =  findobj('type','figure'); % n = length(h); 
%                 for ii = 1:n
%                 name = strcat('GSrejected', num2str(VSDI.ref),'c', num2str(h(ii).Number)); 
%                 saveas(h(ii), fullfile(pathplot, [name, '.jpg']),'jpg')
%                 close(h(ii))
%                 end
%                 clear h n  
%         end
% reject = sort([aberrant_idx0 aberrant_idx1 aberrant_idx2]);
% ######################################################################################

%% Updated: 30/11/20