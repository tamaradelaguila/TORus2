%% TEST CURVAS SOLAPADAS EN NUEVOS FILTROS (from: IGURE FAST SKETCH)

clear

% roi to analyze
roiname1 = 'dm4m_R';
roiname2 = 'dm2_R'; % dm2-R

% ref_wave = 'circ_filt309';
% ref_movie= '_09filt3' ;

% ref_movie= '_17filt5' ;
 ref_movie= '_18filt6' ;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures' ;%@ SET


load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')

% COMMON SETTINGS FOR THE FUNCTION
feedf.window.min = [-100 100]; % 'feed-function' structure
feedf.window.max = [0 1000];
feedf.window.movsum = 50;
feedf.window.basel = [-100 0];
feedf.window.slope=50;
feedf.window.wmean=[0 250];

feedf.noise.fr_abovenoise = 30;
feedf.noise.SDfactor = 2;

% Window For average-based analysis

feedf.window_ave = feedf.window;

feedf.noise_ave = feedf.noise;
feedf.noise_ave.SDfactor = 4;% SET differences


feedf.method = 'movsum';
% END OF SETTINGS


% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');


addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

% end of user_settings

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

reject_on = [1];  %@ SET
% Subsettings:
setting.manual_reject = 0; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET+
setting.force_include = 0; %@ SET



    
    

%% D – CURVAS COND SOLAPADAS POR REGIÓN

% ref_wave = 'circ_filt517';
%     ref_wave = 'filt517';
    ref_wave = 'circ_filt618';

flagcirc = strcmpi(ref_wave(1:4), 'circ');
    

%----------------------------------------------------------------
% @SET: fish + conditions
%----------------------------------------------------------------
for  block = 1%[1:4]% 1:length(fast_condition_list)
    
    % Get selection to be analyzed from structure
    nfish = groupplot{block,1};
    trial_kinds = groupplot{block,3};
    

    [VSDI] = TORus('load',nfish);
    
    if flagcirc
    roi1 = name2idx(roiname1, VSDI.roi.labels_circ);
    roi2 = name2idx(roiname2, VSDI.roi.labels_circ);
    else 
    roi1 = name2idx(roiname1, VSDI.roi.labels);
    roi2 = name2idx(roiname2, VSDI.roi.labels);

    end
    
        selroi = [roi1 roi2];

    
    VSDroiTS = TORus('loadwave',nfish);
    waves = VSDroiTS.(ref_wave).data; %@ SET
    
    
    
    
    
    %----------------------------------------------------------------
    % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
    %----------------------------------------------------------------
    
    rejectidx = [];
    
    if setting.manual_reject
        rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
    end
    
    if setting.GSabsthres_reject
        rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];
        
    end
    
    if setting.GSmethod_reject
        rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];
    end
    
    if setting.force_include
        rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
        
    end
    
    rejectidx = sort(unique(rejectidx));
    
    
    
    
    %----------------------------------------------------------------
    % SELECT CASES
    %----------------------------------------------------------------
    sel_trials= [];
    
    for condi = makeRow(trial_kinds) %to make sure that only the conditions of interest are computed)
        condtrials = makeCol(find(VSDI.condition(:,1)==condi));
        sel_trials  = [sel_trials; condtrials];
    end
    
    sel_trials= sort(sel_trials);
    
        if reject_on
            sel_trials = setdiff(sel_trials, rejectidx);
        end

    % ONSET OF  WAVES FROM CIRCULAR ROIS FROM ONE CONDITION AVERAGE
    
    % build color matrix
    colors = lines(length(trial_kinds));
    
    figure
    ploti = 1;%counter
    
    for nroi = makeRow(selroi)
        
        ax(ploti) = subplot(1,3,ploti); %waves from all conditions in the first column (loop through conditions)
        hold on;
        
        j = 1; %counter for legend labels
        
        for codi = 1:length(trial_kinds)
            code = trial_kinds(codi);
            tricond = intersect(find(VSDI.condition(:,1) == code) , sel_trials);
            codemA(codi)= VSDI.condition(tricond(1),4); %mA corresponding to the code label (for the next subplot)
            roiwave = mean(waves(:,nroi,tricond),3) ;
            
            hold on;
            output = devo_peak2peak(roiwave, VSDI.timebase, feedf.window, feedf.noise, feedf.method, 0, 0);
            
            idx0= dsearchn(VSDI.timebase, 0);
            
            %             waveW = roiwave(idx0:output.peakidx(2));
            %             slopemean(nroi, codi) = mean(diff(waveW));
            
            %             slopemax(nroi,codi) = output.slopemax;
            %             peak(nroi,codi) = output.peakminusbasel;
            
            plot(VSDI.timebase,roiwave, 'color', colors(codi,:), 'linewidth', 3); hold on %  'linewidth', 1.8
            
            %         ylim([-0.2 .3])
            %             ylim(waveslim.ave)
            
            xlim([-300 600])
            
            ylabel('%\Delta F (trials ave)');
            
            title([ VSDI.roi.labels_circ{nroi}])
            
            
            clear output roiwave waveW slopeval
            
            legend_labels{j} =[ num2str(VSDI.condition(tricond(1),4)) 'mA'];
            legend(legend_labels, 'location', 'southeast')
            
            j = j+1;
            
        end %codi
        
        
        hold off
        
        
        ploti = ploti+1;
        
        
    end %  roi
    
    if flagcirc
        ax(3) =  subplot(1,3,3);
        
        roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R, ax(3));
        axis image
        
    else
        ax(3) =  subplot(1,3,3);
        
        roi_preview_multiple(VSDI.crop.preview, VSDI.roi.manual_poly(selroi,:),ax(3) )
        axis image
    end
    %     set(ax, 'ylim', [-0.05 0.32])
    %         set(ax, 'ylim', [-0.05 0.5])
    
    
    
    sgtitle([num2str(VSDI.ref), '(', num2str(trial_kinds(1)), ')' ])
    
    name2save = fullfile(savein,['plot_D_solapadas' num2str(VSDI.ref) '_' num2str(condi(1)) '_' ref_wave 'reject' num2str(reject_on) '.jpg']);
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 8]);
    saveas(gcf,name2save,'jpg')
    close
end
blob()

%% TEST

%% Updated: 12/10/2021
% Last use: