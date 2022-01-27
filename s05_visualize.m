%% s04 VISUALIZE 
clear
user_settings
nfish = 11; %@ SET

VSDI = TORus('load',nfish);

%% 1. WAVEPLOTS OF TRIALS for each roi
% Load data
VSDroiTS = TORus('loadwave',nfish);
waves = VSDroiTS.filt306.data; %@ SET

% Make waveplots
titulo = strcat('test_fish',num2str(VSDroiTS.ref)); %@ SET
close all
plot_wavesplot(waves, VSDroiTS.roi.labels, titulo, VSDI.nonanidx, 20);

% Save and close
path2save = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot' ; %@ SET
name2save = strcat('waveplot', num2str(VSDroiTS.ref));  %@ SET (if you want)

save_currentfig(path2save, name2save)
blob()

%% 2. ERPs

VSDroiTS = TORus('loadwave',nfish);

%@ SET: what to plot:
timeseries= VSDroiTS.filt306.data;

roi2plotidx = 1:14;%[1 3 5 7 11];
trials2plot = VSDI.nonanidx; %@ SET 
% trials2plot = choosetrials(VSDI.condition, [1 2]);
% plot_ERPs(timeseries, VSDI.timebase,trials2plot, roi2plotidx , VSDroiTS.roi.labels, titletxt, stim_dur); 

% Automatic extraction and saving
      code2loop = sort(unique(VSDI.condition(:,1)));
      code2loop = code2loop(~isnan(code2loop));
    
        for ii =1:length(code2loop) %loop through all condition codes
        code = code2loop(ii);
        trials2plot = find(VSDI.condition(:,1) == code); %all trials will have the same charact
        tcode = num2str(code);
        cond_def = strcat('code#',tcode(1),'-',tcode(2),'-',num2str(VSDI.list(trials2plot(1)).mA),'mA');
        titletxt = strcat(num2str(VSDroiTS.ref),cond_def);

        stim_dur = VSDI.list(trials2plot(1)).Sdur; % take the duration of the stimulus from one of the selected trials
        yaxis = [-0.2 0.7];
        plot_ERPs(timeseries, VSDI.timebase,trials2plot, roi2plotidx , VSDroiTS.roi.labels, titletxt, stim_dur, yaxis); 

            % Save and close
            path2save = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps' ; %@ SET
            name2save = strcat(num2str(VSDroiTS.ref),'pp09_erps_',cond_def,'circ') ;  %@ SET (if you want)
            ...name2save = strcat(num2str(VSDroiTS.ref),'erps_',cond_def) ;  %@ SET (if you want)

            save_currentfig(path2save, name2save)

        end
        
        

%% 2.extra ERPs + preview region

clearvars
user_settings

 nfish = 11; 
    VSDI = TORus('load',nfish);
    VSDroiTS = TORus('loadwave',nfish);
    
    ....@ SET: what to plot:
%     timeseries= VSDroiTS.circ_filt309.data;
    timeseries= VSDroiTS.circ_filt306.data;

%     timeseries= VSDroiTS.filt306.data;
        
    roi2plotidx =[1:14];%[1 3 5 7 11];
    % trials2plot = choosetrials(VSDI.condition, [1 2]);
    % plot_ERPs(timeseries, VSDI.timebase,trials2plot, roi2plotidx , VSDroiTS.roi.labels, titletxt, stim_dur);
    
    ....Automatic extraction and saving
    %       code2loop = {'200', '201', '202' , '203', '300', '301', '302','303'}; % @SET fish 23
    code2loop = sort(unique(VSDI.condition(:,1)));
    code2loop = code2loop(~isnan(code2loop));
    
    for ii = makeRow(code2loop) %loop through all condition codes
        code = num2str(ii);
        trials2plot = find(VSDI.condition(:,1) == ii); %all trials will have the same charact
        %         cond_def =
        %         strcat('code#',code(1),'-',code(2),'-',code(3),'(',num2str(VSDI.list(trials2plot(1)).mA),'mA)');
        %         %only for 3-letters code
        cond_def = strcat('code#',code,'(',num2str(VSDI.list(trials2plot(1)).mA),'mA)');
        
        titletxt = strcat(num2str(VSDroiTS.ref),cond_def);
        
        stim_dur = VSDI.list(trials2plot(1)).Sdur; % take the duration of the stimulus from one of the selected trials
        yaxis = [-0.1 0.5];
        plot_ERPs_plusroi(timeseries, VSDI.timebase,trials2plot, roi2plotidx , VSDroiTS.roi.labels, titletxt, stim_dur, VSDI.crop.preview, VSDI.roi.manual_poly, yaxis);
        
        
        ... Save and close
        path2save = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps' ; %@ SET
        name2save = strcat(num2str(VSDroiTS.ref),'pp06_erps_',cond_def,'circ') ;  %@ SET (if you want)

        %name2save = strcat(num2str(VSDroiTS.ref),'erps_',cond_def) ;  %@ SET (if you want)
        
        save_currentfig(path2save, name2save)
        
    end
    clearvars -except path nfish
    %         end %for nfish = [5]
    
    %% 3. FRAME from movie
    
    %@ SET...
    % ... what and how to plot
    movie_ref = '_06filt3'; % input movie
    frame2plot = 66; %idx
    trial2plot = 2; %idx
    act_clim= [-0.2 0.2]; %coloraxis of the shown colors
    % ...noise threshold settings
    baseline = 1:20; %range of frames to consider baseline (idx)
    SDfactor = 4;
    %...absolute threshold settings
    method = 'absvalue';
    value = 0.07;
    % Load movie
    VSDmov = TORus('loadmovie',nfish,movie_ref);
    movie = VSDmov.data(:,:,1:end-1,trial2plot); %movie to plot
    
    %@ SET which threshold to use (and leave the other unmuted)
    % NOISE THRESHOLD:
    % [movie_thres, alphachan, ~] = movie_noisethresh(movie,baseline,SDfactor, 0);
    
    % ABSOLUTE THRESHOLD
    [movie_thres, alphachan] = movie_absthresh(movie,method,value);
    % NOTE_DEV: comprobar c�mo tiene que ser alphachannel
    
    % OVERLAID PLOT
    background = VSDI.backgr(:,:,trial2plot);
    % frame = movie(:,:,frame2plot);
    frame = movie_thres(:,:,frame2plot);
    framealpha = alphachan(:,:,frame2plot);
    
    ax1 = subplot(1,1,1);
    title(['']);
    plot_framesoverlaid(frame,background, ~framealpha ,0, ax1, act_clim, 1);
    % imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar
    ax1.XTick = []; ax1.YTick = [];
    
    % DEV NOTE (still to do):  imlpement tiles function with the overimposed
    % colors
    
    %% 4 - TILES
    movie_ref = '_06filt3'; % input movie
    VSDmov = TORus('loadmovie',nfish,movie_ref);
    
    tileset.start_ms = -100; % time in ms for first tile
    tileset.end_ms = 1000;
    tileset.nrows = 6;
    tileset.ncols = 8;
    %           tileset.climsat = 0.8 ; %colormap limit would be the 80% of the max/min value
    tileset.clims = [-0.2 0.2];
    tileset.thresh = [-0.05 0.05];
    tileset.time2plot = 0; %select time (ms)
    tileset.x = 35;
    tileset.y = 35;
    
    movie2plot = VSDmov.data(:,:,:,2);
    plot_tilemovie(movie2plot, VSDI.timebase, tileset);
    
    %% 5 - OVERLAID 12 TILES : from a trial
    movie_ref = '_06filt3'; % input movie
    VSDmov = TORus('loadmovie',nfish,movie_ref);
    trial = 30;
    
    tileset.start_ms = -100; % time in ms for first tile
    tileset.end_ms = 800;
    tileset.clims = [-0.9 0.9];
    tileset.thresh = [-0.1 0.1];
    movie2plot = VSDmov.data(:,:,:,trial);
    
    plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
    
    %% 6 - OVERLAID 12 TILES : from an average
    
    movie_ref = '_06filt3'; % input movie
    VSDmov = TORus('loadmovie',nfish,movie_ref);
    
    [idxA] = find(VSDI.condition(:,1)==102);
    [idxB] = find(VSDI.condition(:,1)==100);
    
    
    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,idxA),4) ;
    
    movie2plot_dif = mean(VSDmov.data(:,:,:,idxA),4) - mean(VSDmov.data(:,:,:,idxB),4) ;
    movie2plot_dif(:,:,end) = movie2plot(:,:,end); %substitute the substracted background with a normal one
    
    % settings
    tileset.start_ms = -100; % time in ms for first tile
    tileset.end_ms = 800;
    %           tileset.clims = [-0.9 0.9];
    tileset.clims = [-0.3 0.3];
    tileset.thresh = [-0.1 0.1];
    
    plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
    titulo = ['AVERAGE cond.101', '(',num2str(VSDI.ref) ,')'];
    sgtitle(titulo)
    
    
    plot_tilemovie12frames(movie2plot_dif, VSDI.timebase, tileset);
    titulo2 = ['AVERAGE cond.101 - control', '(',num2str(VSDI.ref) ,')'];
    sgtitle(titulo2)
    
    
    %% 6 - OVERLAID ANY NUMBER OF TILES : from a movie
    
    movie_ref = '_12filt5'; % input movie
    VSDmov = TORus('loadmovie',nfish,movie_ref);
    
    [idxA] = find(VSDI.condition(:,1)==404);
    [idxB] = find(VSDI.condition(:,1)==400);
    
    
    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,idxA),4) ;
    
    movie2plot_dif = mean(VSDmov.data(:,:,:,idxA),4) - mean(VSDmov.data(:,:,:,idxB),4) ;
    movie2plot_dif(:,:,end) = movie2plot(:,:,end); %substitute the substracted background with a normal one
    
    % settings
    tileset.start_ms = -100; % time in ms for first tile
    tileset.end_ms = 800;
    %           tileset.clims = [-0.9 0.9];
    tileset.clims = [-0.3 0.3];
    tileset.thresh = [-0.1 0.1];
        

    plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
    plot_tilemovie12frames_devo(movie2plot, VSDI.timebase, tileset);
    
roi2draw = VSDI.roi.manual_poly([1,2,7,8,13,14]);
    plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset, roi2draw);
    plot_tilemovie_custom(movie2plot, VSDI.timebase, tileset);

    titulo = ['AVERAGE cond.101', '(',num2str(VSDI.ref) ,')'];
    sgtitle(titulo)
    
    
    plot_tilemovie12frames(movie2plot_dif, VSDI.timebase, tileset);
    titulo2 = ['AVERAGE cond.101 - control', '(',num2str(VSDI.ref) ,')'];
    sgtitle(titulo2)
    
    
    %%
    clear
    user_settings
    
    
    for nfish = [10]
        VSDI = TORus('load',nfish);
        movie_ref = '_06filt3'; % input movie
        VSDmov = TORus('loadmovie',nfish,movie_ref);
        
        %@ SET: what to plot:
        code2loop = {'2101', '2100', '3101' , '3100'}; % @SET fish 23
        %       code2loop = {'100', '101', '102' , '103'}; % @SET
        
        for ii =1:length(code2loop) %loop through all condition codes
            code = code2loop{ii};
            
            
            
            [idxA] = find(VSDI.condition(:,1)==101);
            [idxB] = find(VSDI.condition(:,1)==contr);
            
            %to plot single trial
            movie2plot = mean(VSDmov.data(:,:,:,idxA),4) ;
            
            movie2plot_dif = mean(VSDmov.data(:,:,:,idxA),4) - mean(VSDmov.data(:,:,:,idxB),4) ;
            movie2plot_dif(:,:,end) = movie2plot(:,:,end); %substitute the substracted background with a normal one
            
            % settings
            tileset.start_ms = -100; % time in ms for first tile
            tileset.end_ms = 800;
            %           tileset.clims = [-0.9 0.9];
            tileset.clims = [-0.3 0.3];
            tileset.thresh = [-0.1 0.1];
            
            plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
            titulo = ['AVERAGE cond.101', '(',num2str(VSDI.ref) ,')'];
            sgtitle(titulo)
            
            
            plot_tilemovie12frames(movie2plot_dif, VSDI.timebase, tileset);
            titulo2 = ['AVERAGE cond.101 - control', '(',num2str(VSDI.ref) ,')'];
            sgtitle(titulo2)