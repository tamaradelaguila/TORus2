%% PLOT WAVES AND STANDART DEVIATION
clear

% ...........................................
% DOLOR
% ...........................................

% Dm4
% ...........................................
dm4 = load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/meanwave_std/group19_dolorn4_normal_dm4.csv');
dm4 = dm4';

timeb4 = dm4(1,:);
dm4 = dm4(2:end, :);

SEM4 = std(dm4)/sqrt(size(dm4,1));               % Standard Error

% PLOT STANDART ERROR
figure
lower4 =mean(dm4)  - std(dm4); 
upper4 = mean(dm4) + std(dm4); 
ciplot(lower4,upper4,timeb4,'b'); hold on
plot(timeb4, mean(dm4), 'b', 'linewidth', 3'); hold on

% Dm2
% ...........................................

dm2 = load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/meanwave_std/group19_dolorn4_normal_dm2.csv');
dm2 = dm2';

timeb2 = dm2(1,:);
dm2 = dm2(2:end, :);

SEM2 = std(dm4)/sqrt(size(dm2,1));               % Standard Error

% PLOT STANDART ERROR
lower2 =mean(dm2)  - std(dm2); 
upper2 = mean(dm2) + std(dm2); 
ciplot(lower2,upper2,timeb2,'r'); hold on
plot(timeb2, mean(dm2), 'r', 'linewidth', 3')
title('dolor_n4')
xlabel('ms')
ylabel('mean+-std')

% 
% 
% % PLOT SEM
% lower =mean(data) + SEM; 
% ciplot(lower,upper,x,colour);


% ...........................................
% CHORRITO
% ...........................................
clear

% Dm4
% ...........................................
dm4 = load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/meanwave_std/group16_chorriton3_normal_dm4.csv');
dm4 = dm4';

timeb4 = dm4(1,:);
dm4 = dm4(2:end, :);

SEM4 = std(dm4)/sqrt(size(dm4,1));               % Standard Error

% PLOT STANDART ERROR
figure
lower4 =mean(dm4)  - std(dm4); 
upper4 = mean(dm4) + std(dm4); 
ciplot(lower4,upper4,timeb4,'b'); hold on
plot(timeb4, mean(dm4), 'b', 'linewidth', 3'); hold on

% Dm2
% ...........................................

dm2 = load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/meanwave_std/group16_chorriton3_normal_dm2.csv');
dm2 = dm2';

timeb2 = dm2(1,:);
dm2 = dm2(2:end, :);

SEM2 = std(dm4)/sqrt(size(dm2,1));               % Standard Error

% PLOT STANDART ERROR
lower2 =mean(dm2)  - std(dm2); 
upper2 = mean(dm2) + std(dm2); 
ciplot(lower2,upper2,timeb2,'r'); hold on
plot(timeb2, mean(dm2), 'r', 'linewidth', 3')
title('dolor_n4')
xlabel('ms')
ylabel('mean+-std')


%% CONFIDENCE INTERVAL
data = load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/05_figure_definit/meanwave_std/group19_dolorn4_normal_dm4.csv');
data = data';

for wavi = 1:size(data,1)
wave = data(wavi, :); 
SEM = std(wave)/sqrt(length(wave));               % Standard Error
ts = tinv([0.025  0.975],length(wave)-1);      % T-Score
CI = mean(wave) + ts*SEM;    

end


lower =mean(wave) + SEM; 
ciplot(lower,upper,x,colour);

%% Created: 05/10/22
