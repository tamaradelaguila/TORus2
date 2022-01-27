% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% BLANK SUBSTRACTION OF THE WAVES

clear
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')

% ///////////////////////////////////////////////////////////
% SETTINGS

 sourcedata = 'normal';
% sourcedata = '%F';
% sourcedata = 'blank-s'; % BLANK-SUBSTRACTION
% sourcedata = '%deltaF blank-s';

% selroinames = {'dm4m_R',  'dm2_R'};
selroinames = {'dm4m_R',  'dm2_R', 'dm3_R' ,'dm1_R','dldm_R'};%dm3,

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_17filt5' ;
% ref_movie= '_12filt5' ;

activ_unit = 'diffF'; % @ MANUALLY SET

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures' ;%@ SET

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis


% ///////////////////////////////////////////////////////////////

% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');

addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));
% END USER_SETTINGS

for reject_on = [1]  %@ SET
    
    %----------------------------------------------------------------
    % @SET: REJECT SETTINGS
    %----------------------------------------------------------------
    
    % Subsettings:
    setting.manual_reject = 1; %@ SET
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET+
    setting.force_include = 0; %@ SET
    
    
    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    
    for suji  =  1:size(groupplot,1)
        
        nfish = groupplot{suji,1};
        cond_lohi = groupplot{suji,3};
        
        VSDI = TORus('load', nfish);
        VSDmov = TORus('loadmovie',nfish,ref_movie);
        movies = VSDmov.data ;
        F0 = VSDmov.F0;
        
        %----------------------------------------------------------------
        % CONTROL ROI PICTURE
        %----------------------------------------------------------------
        %     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R);
        %     title([num2str(VSDI.ref) 'roi preview:' selroinames{:}])
        %     saveas(gcf, [num2str(VSDI.ref)'roipreview'], 'jpg')
        
        %----------------------------------------------------------------
        % SELECT ROI
        %----------------------------------------------------------------
        
        switch roikind
            case 'circle'
                selroi =name2idx(selroinames, VSDI.roi.labels_circ);
                roilabels = VSDI.roi.labels_circ;
            case 'anat'
                selroi =name2idx(selroinames, VSDI.roi.labels);
                roilabels = VSDI.roi.labels;
        end
        
        %----------------------------------------------------------------
        % GET INDEXES OF TIMERANGE
        %----------------------------------------------------------------
        idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
        idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values
        
        timebase_adj = VSDI.timebase(idxrange);
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
        
        % -------------------------------------------
        % FIRST LOOP THROUGH CONDITIONS AND ROIS TO GET MEASURES AND LATENCY RESPECT TO THE MAX PEAK (of all regions, i.e., a
        % common threshold to all)
        % -------------------------------------------
        % we get and store the value in this first loop to get the max
        % value, and in the second loop we use that max-val as
        % threshold
        
        
        for condition =  makeRow(cond_lohi)
            
            
            cond_blank = force0ending(condition);
            % -------------------------------------------
            % SELECT CASES  AND AVERAGE MOVIE
            % -------------------------------------------
            sel_trials = find(VSDI.condition(:,1)==condition);
            sel_blank = find(VSDI.condition(:,1)==cond_blank);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
                sel_blank = setdiff(sel_blank, rejectidx);
                
            end
            
            for roi_i =  makeRow(selroi)
                
                % STORE THE WAVES TO PLOT
                
                switch roikind
                    case 'circle'
                        roimask = VSDI.roi.circle.mask(:,:,roi_i);
                        
                        %             anamask = VSDI.roi.manual_mask(:,:,1);
                    case 'anat'
                        roimask = VSDI.roi.manual_mask(:,:,1);
                end
                
                % --------------------------------------------------------------------------
                % STORE WAVES. Use timerange set
                
                
                switch sourcedata
                    
                    case 'normal'
                        movieave = mean(movies(:,:,idxrange,sel_trials),4);
                        movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                        ...(that normally corresponds to the background)
                            
                    
                    roiwave =  roi_TSave(movieave,roimask);
                    refcase = '';
                    maxval = max(roiwave(:));

                    
                    case  '%F'
                        
                        % Step-1: calculate %F trial-wise
                        for triali =1:size(movies,4)
                            movie = movies(:,:,idxrange,triali);

                            movie(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                            trialwave(:,triali) =   roi_TSave_percF_roiwise(movie,roimask, F0(:,:,triali))*1000; % scaling factor
                        end
                        % Step-2: average waves from selected trials
                        roiwave  = mean(trialwave(:,sel_trials),2);
                        refcase = '%F';
                        
                        % Step-3 : get non %F wave to get the color limit
                        % for the tile
                        movieave = mean(movie, 4); 
                        localwave =  roi_TSave(movieave,roimask);
                        maxval = max(localwave(:)); 


                    case  'blank-s'
                        movieave = mean(movies(:,:,idxrange,sel_trials),4);
                        movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                        
                        movieblank = mean(movies(:,:,idxrange,sel_blank),4);
                        movieblank(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                        blankwave = roi_TSave(movieblank,roimask);
                        roiwave =  roi_TSave(movieave,roimask)-blankwave;
                        refcase = 'blankS';
                        %             case '%deltaF blank-s'
                        %                 refcase = '%deltaF blanks';
                        maxval = max(roiwave(:));
                        
                end %case
                roiwaves(:,roi_i) = roiwave;
                

            end % roi
            
            %% PLOT CONDITION-WISE
            

            % PLOT BRAINPREVIEW
            ax1 = subplot(2,3,3);

            switch roikind
                case 'circle'
                    % selroi =name2idx(selroinames, VSDI.roi.labels_circ);
                    % roilabels = VSDI.roi.labels_circ;
                    centers = VSDI.roi.circle.center(selroi, :) ;
                    
                    roicirc_preview_multiple(VSDI.crop.preview, centers, VSDI.roi.circle.R, ax1);
                    
                case 'anat'
                    % selroi =name2idx(selroinames, VSDI.roi.labels);
                    % roilabels = VSDI.roi.labels;
                    roi_preview_multiple(VSDI.crop.preview, VSDI.roi.manual_poly(selroi,:), ax1);
                    
            end
            
            
            % PLOT WAVES
            subplot(2, 3, [1 2 4 5])
            waves2plot = roiwaves(:,selroi);
            plot(timebase_adj, waves2plot, 'linewidth', 2);
            legend(selroinames)
            title(num2str(suji))
            
            localref = [num2str(VSDI.ref) '- cond=' num2str(condition) '-' ref_movie '-' refcase  '-reject' num2str(reject_on) ]; 
            
            switch sourcedata
                case 'normal'
                 ylabel(activ_unit)

                    
                case  '%F'
                    
                   ylabel(['%F (from' activ_unit ')'])

                case  'blank-s'
                  ylabel([activ_unit '(blank substracted)'])

                    
                    
            end %case

            
            sgtitle(localref)
            
            savename = ['TILEWAVES' localref   '-' num2str(numel(selroi)) roikind 'ROI - PART2 WAVES'];

            saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
            close

            
          % PLOT TILES
          movieplot = mean(movies(:,:,idxrange,sel_trials),4);
          
          
          tileset.start_ms = 0; % time in ms for first tile       
          tileset.end_ms = 800;
          tileset.time2plot = 0; %select time (ms)
          tileset.x = 35; 
          tileset.y = 35; 
          tileset.clims = [-maxval maxval] ;
          tileset.thresh = [-maxval/6 maxval/6];
          tileset.nrowcol = [6 4];
          tileset.backgr = VSDI.crop.preview;

          plot_tilemovie_custom(movieplot, timebase_adj, tileset,[])
           
          
            
            localref = [num2str(VSDI.ref) '- cond=' num2str(condition) '-' ref_movie '-' refcase  '-reject' num2str(reject_on) ]; 
            sgtitle(localref)
            
            savename = ['TILEWAVES' localref   '-' num2str(numel(selroi)) 'ROI - PART1 TILES'];

            saveas(gcf, fullfile(savein, [savename '.jpg']), 'jpg')
            close
        end % condi
        
        blob()
    end % for suji
end % for reject_on