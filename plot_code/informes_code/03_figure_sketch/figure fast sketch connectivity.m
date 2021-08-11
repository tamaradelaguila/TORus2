%% FUNCTIONAL CONNECTIVITY TEST

% OPTION 1 - Get average ROI timeseries from each trial, calculate the
% correlation matrix and then average those.

% OPTION 2 - Get average movie, extract the ROI timeserie and calculate the
% correlation matrix

% OPTION 3 -Get timeseries from each movie, average them and get the
% connectivity matrix

%% OPTION 1 - average correlation matrix from each trial
clear
ref_wave = 'circ_filt309';
ref_movie= '_09filt3' ;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')


% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

% end of user_settings

% selroinames = {'dm4m_R', 'dm2_R', 'dm3_R', 'dldm_R', 'dm1_R'}; % @SET
selroinames = {'dm4_R', 'dm4_L','dm2_R', 'dm2_L','dm3_R', 'dm3_L','dldm_R',  'dldm_L','dm1_R',  'dm1_L'}; ; % @SET


for block = 1:length(fast_condition_list)
    
    nfish = fast_condition_list{block,1};
    
    trial_kinds = fast_condition_list{block,2};
    cond_def =fast_condition_list{block,3};
    
    VSDI = TORus('load', nfish);
    VSDroiTS = TORus('loadwave',nfish);
    
    selroi =name2idx(selroinames, VSDI.roi.labels_circ); %idx of selected subset
    waves = VSDroiTS.(ref_wave).data(101:151,selroi,:); % waves is already a subset 

    nplot = length(trial_kinds);
    nroi = size(waves,2);

    %--------------------------------------
    % 1. GET CONNECTIVITY MATRIX FOR EACH CONDITION
    %--------------------------------------
    j=1;
    for condi = makeRow(trial_kinds)
        [idxA] = find(VSDI.condition(:,1)==condi);
        
        for triali = makeRow(idxA)
            
            % APPLY FUNCTION
            tempFCmatrix(:,:,triali) =FConn_matrix(waves(:,:,triali));
            
            % ADJUSTMENTS
        end % triali
        
        conM(:,:,j) = mean(tempFCmatrix, 3);
        
        mA(j)=VSDI.condition(idxA(1),4);
        
        j = j+1;
        
    end % condi
    
    %--------------------------------------
    % 2. PLOT CONNECTIVITY MATRIX
    %--------------------------------------
    if nplot <= 4
       n = [2 2];
    elseif nplot >4
       n = [2 3] ;
    end
    
    for ploti = 1:size(conM,3)
        
        ax(ploti) = subplot(n(1),n(2),ploti);
        imagesc(conM(:,:,ploti));
        axis square
        colormap(parula);     colorbar
        
        xticks(1:nroi)
        xticklabels (selroinames); xtickangle(90)
        yticks(1:nroi)
        yticklabels (selroinames)
        title(['c=' num2str(trial_kinds(ploti)) '(' num2str(mA(ploti))  'mA)'])
    end
%             set(ax, 'clim' ,[-1 1], 'colormap', flipud(RdBu))
            ccolor = parula; ccolor(1,:) = [.3 .3 .3];
            set(ax, 'clim' ,[0 .8], 'colormap', ccolor)

    sgtitle ([num2str(VSDI.ref),': connect mat (\rho )'])
    
    %--------------------------------------
    % 3. SAVE IMAGE
    %--------------------------------------
    name2save = fullfile(savein,['plot_F_connM_opt1',num2str(VSDI.ref),'cod' , '.jpg']);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 10]);
    saveas(gcf,name2save,'jpg')
    close
    
    
end %block

blob()

%% OPTION 2 - correlation matrix from averaged movie
clear
ref_wave = 'circ_filt309';
ref_movie= '_09filt3' ;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_plot' ;%@ SET
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/fast_condition_list.mat')


% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

% end of user_settings

% selroinames = {'dm4_R','dm4m_R', 'dm2_R', 'dm3_R', 'dldm_R', 'dm1_R'}; % @SET
selroinames = {'dm4_R', 'dm2_R', 'dm3_R', 'dldm_R', 'dm1_R', 'dm4_L', 'dm2_L', 'dm3_L', 'dldm_L', 'dm1_L'}; ; % @SET


for block = 1:length(fast_condition_list)
    
    nfish = fast_condition_list{block,1};
    
    trial_kinds = fast_condition_list{block,2};
    cond_def =fast_condition_list{block,3};
    
    VSDI = TORus('load', nfish);
    VSDmov = TORus('loadmovie',nfish,ref_movie);
    
    selroi =name2idx(selroinames, VSDI.roi.labels_circ);

    nplot = length(trial_kinds);
    nroi = size(waves,2);
    
    %--------------------------------------
    % 1. GET AVERAGE MOVIE FOR EACH CONDITION
    %--------------------------------------
    j=1;
    for condi = makeRow(trial_kinds)
        [idxA] = find(VSDI.condition(:,1)==condi);
        
        for triali = makeRow(idxA)
            
            
            % ADJUSTMENTS
        end % triali
        
        conM(:,:,j) = mean(tempFCmatrix, 3);
        
        mA(j)=VSDI.condition(idxA(1),4);
        
        j = j+1;
        
    end % condi
    
    %--------------------------------------
    % 2. GET CONNECTIVITY MATRIX
    %--------------------------------------
            % APPLY FUNCTION
            tempFCmatrix(:,:,triali) =FConn_matrix(waves(:,selroi,triali));

    
    %--------------------------------------
    % 3 PLOT CONNECTIVITY MATRIX
    %--------------------------------------
    if nplot <= 4
       n = [2 2];
    elseif nplot >4
       n = [2 3] ;
    end
    
    for ploti = 1:size(conM,3)
        
        ax(ploti) = subplot(n(1),n(2),ploti);
        imagesc(conM(:,:,ploti));
        axis square
        colormap(parula);     colorbar
        
        xticks(1:nroi)
        xticklabels (selroinames); xtickangle(90)
        yticks(1:nroi)
        yticklabels (selroinames)
        title(['c=' num2str(trial_kinds(ploti)) '(' num2str(mA(ploti))  'mA)'])
    end
%             set(ax, 'clim' ,[-1 1], 'colormap', flipud(RdBu))
            ccolor = parula; ccolor(1,:) = [.3 .3 .3];
            set(ax, 'clim' ,[0 .8], 'colormap', ccolor)

    sgtitle ([num2str(VSDI.ref),': connect mat (\rho )'])
    
    %--------------------------------------
    % 4 SAVE IMAGE
    %--------------------------------------
    name2save = fullfile(savein,['plot_F_connM_opt2',num2str(VSDI.ref),'cod' , '.jpg']);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 10 10]);
    saveas(gcf,name2save,'jpg')
    close
    
    
end %block

blob()

