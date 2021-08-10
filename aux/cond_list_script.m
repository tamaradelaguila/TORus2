%% condition list to loop through to get different conditions from each fish and block

%% FISH 210320
n= 1;
condition_list{n,1} =1 ; %FISH 
condition_list{n,2} = 110:113; 
condition_list{n,3} = 'intrad 03ms';

n = n+1; 
condition_list{n,1} =1 ; %FISH 
condition_list{n,2} = 300:303; 
condition_list{n,3} = 'tren largo p10ms';

%% FISH 210323
n = n+1; 
condition_list{n,1} = 2; %FISH 
condition_list{n,2} = 200:203; 
condition_list{n,3} = 'tren p03ms';

n = n+1; 
condition_list{n,1} = 2; %FISH 
condition_list{n,2} = 300:303; 
condition_list{n,3} = 'tren largo p10ms';

%% FISH 210325
% n = n+1; 
% condition_list{n,1} = 3; %FISH 
% condition_list{n,2} = 100:103; 
% condition_list{n,3} = 'pulso 03ms';

n = n+1; 
condition_list{n,1} = 3; %FISH 
condition_list{n,2} = 300:303; 
condition_list{n,3} = 'tren largo p10ms';

%% FISH 210409

n = n+1; 
condition_list{n,1} = 4; %FISH 
condition_list{n,2} = 100:103; 
condition_list{n,3} = 'pulso 03ms';

%% FISH 210411

% n = n+1; 
% condition_list{n,1} = 5; %FISH 
% condition_list{n,2} = 1101; 
% condition_list{n,3} = 'pulso 03ms no aleat';

n = n+1; 
condition_list{n,1} = 5; %FISH 
condition_list{n,2} = 2100:2101; 
condition_list{n,3} = 'pulso 03ms +control';

n = n+1; 
condition_list{n,1} = 5; %FISH 
condition_list{n,2} = 3100:3101; 
condition_list{n,3} = 'pulso 03ms +control post sobredosis';

%% FISH 210412

% n = n+1; 
% condition_list{n,1} = 6; %FISH 
% condition_list{n,2} = [1100 1101 1401]; 
% condition_list{n,3} = 'control pulso015 tren';
% 
% n = n+1; 
% condition_list{n,1} = 6; %FISH 
% condition_list{n,2} = [2000:2003]; 
% condition_list{n,3} = 'control 3 pulsitos 03ms';

n = n+1; 
condition_list{n,1} = 6; %FISH 
condition_list{n,2} = [3000:3002]; 
condition_list{n,3} = 'control,tono,tren premorfina';

n = n+1; 
condition_list{n,1} = 6; %FISH 
condition_list{n,2} = [4000:4002]; 
condition_list{n,3} = 'control,tono,tren postmorfina';


%% FISH 210430

n = n+1; 
condition_list{n,1} = 8; %FISH 
condition_list{n,2} = [1000:1003]; 
condition_list{n,3} = 'tono pulso015';

n = n+1; 
condition_list{n,1} = 8; %FISH 
condition_list{n,2} = [2000:2003]; 
condition_list{n,3} = 'pulso015';

%% FISH 210508

n = n+1; 
condition_list{n,1} = 9; %FISH 
condition_list{n,2} = [1000:1003]; 
condition_list{n,3} = 'tono pulso015';

%% FISH 210509

n = n+1; 
condition_list{n,1} = 10; %FISH 
condition_list{n,2} = [1000:1003]; 
condition_list{n,3} = 'tono pulso015';

n = n+1; 
condition_list{n,1} = 10; %FISH 
condition_list{n,2} = [2000:2003]; 
condition_list{n,3} = 'tono pulso015';

%% FISH 210521

n = n+1; 
condition_list{n,1} = 11; %FISH 
condition_list{n,2} = [400:403]; 
condition_list{n,3} = 'pulso015';

%% FISH 210522

n = n+1; 
condition_list{n,1} = 12; %FISH 
condition_list{n,2} = [400:403]; 
condition_list{n,3} = 'pulso015';

save('/home/tamara/Documents/MATLAB/VSDI/TORus/data/condition_list.mat', 'condition_list')