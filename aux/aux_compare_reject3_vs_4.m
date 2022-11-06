ref = 'filt6_18';
ref2 = 'filt6_21' ;
nfish = 1; 
VSDI = TORus('load', nfish)


pathmat = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataVSDI/';
load(fullfile(pathmat, 'TORus_GSwaves'))

% -------------------------------
% plot individual waves to see
% -------------------------------
for n = makeRow(VSDI.nonanidx)
close
% subplot(1,3,1)
% plot(GS.(ref).pF{nfish}(:,n)*100, 'Displayname', ref)
% title([ref])
% 
% subplot(1,3,2)
plot(GS.(ref2).wave{nfish}(:,n), 'Displayname', ref2)
ylim([-.1 0.35])
title(ref2)

% subplot(1,3,3)
% plot(GS.(ref).pF{nfish}(:,n)*100, 'Displayname', [ref '*100'])
% hold on 
% plot(GS.(ref2).wave{nfish}(:,n), 'Displayname', ref2)
% ylim([0 0.35])
% legend('location', 'southeast')
% title(num2str(n))

pause
end


% -------------------------------
% COMPARE reject3 VS reject4
% -------------------------------
nfish =12; 
VSDI = TORus('load', nfish);

% GSabs
m1 = ['reject3-GSabs = ' num2str( numel(VSDI.reject3.GSabs025))];
m2 = ['reject4-GSabs = ' num2str( numel(VSDI.reject4.GSabs025))];
m3 = ['GSabs overlap in:' num2str(numel(intersect(VSDI.reject3.GSabs025, VSDI.reject4.GSabs025)))];
    m4 =  setdiff(VSDI.reject4.GSabs025, VSDI.reject4.GSabs025);

    m5 =  setdiff(VSDI.reject3.GSabs025, VSDI.reject4.GSabs025);

disp(VSDI.ref)
disp(m1)
disp(m2)
disp(m3)


disp('indexes in reject3 not in reject4:')
disp(m5)

disp('indexes in reject4 not in reject3:')
disp(m4)

%GSdeviation
m1 = ['reject3-GSdeviat2sd = ' num2str( numel(VSDI.reject3.GSdeviat2sd))];
m2 = ['reject4-GSdeviat2sd = ' num2str( numel(VSDI.reject4.GSdeviat2sd))];
m3 = ['GSdeviat2sd overlap in:' num2str(numel(intersect(VSDI.reject3.GSdeviat2sd, VSDI.reject4.GSdeviat2sd)))];
    m4 =  setdiff(VSDI.reject4.GSdeviat2sd, VSDI.reject3.GSdeviat2sd);

    m5 =  setdiff(VSDI.reject3.GSdeviat2sd, VSDI.reject4.GSdeviat2sd);

disp(VSDI.ref)
disp(m1)
disp(m2)
disp(m3)


disp('indexes in reject3 not in reject4:')
disp(m5)

disp('indexes in reject4 not in reject3:')
disp(m4)

% -------------------------------
% PLOT INDIVIDUAL WAVE TO CHECK
% -------------------------------
idx = 206;
waves = GS.filt6_21.wave{nfish};
plot(waves,'color',[.5 .5 .5]); hold on
plot(waves(:,idx), 'r', 'linewidth', 1.3)
ylim([-.1 0.35])
