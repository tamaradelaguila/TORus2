% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% CALCULATES THE MAPS TO LATER PLOT THEM IN THE SAME SCALE TO DETECT
% OUTLIERS 


clear; W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W); clear W

% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new.mat') % ATTENTION: select cases manually (there are repetitions)

% ///////////////////////////////////////////////////////////

    %----------------------------------------------------------------
    % @SET: SELECT FISH AND CONDITIONS
    %----------------------------------------------------------------
    
    
    nfish = 10;
    cond_codes = [2003];
    
ref_movie= '_18filt6';% '_17filt5' ; '_18filt6'
% ref_movie= '_12filt5' ;

reject_on =0;  %@ SET


% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 600]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-25 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350];

% feedf.noise.fr_abovenoise = 30;
% feedf.noise.SDfactor = 2;
% feedf.noise.basel = [-200 0]; %it won't be used by the function, but will be used to manually get the baseline
feedf.method = 'movsum';

% flagslope = 1;

% Params
slope.window = [0, 200]; %ms


% ///////////////////////////////////////////////////////////////

% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');

addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));

    
    %----------------------------------------------------------------
    % @SET: REJECT SETTINGS
    %----------------------------------------------------------------
    
    % Subsettings:
    setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET+
    setting.force_include = 0; %@ SET
    

    %----------------------------------------------------------------
    % LOAD DATA
    %----------------------------------------------------------------
    
   
        
        VSDI = TORus('load', nfish);
        VSDmov = TORus('loadmovie',nfish,ref_movie);

        %----------------------------------------------------------------
        % GET INDEXES OF TIMERANGE
        %----------------------------------------------------------------
        idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
        idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values
        
        timebase_adj = VSDI.timebase(idxrange);
        
        
        % get indexes for slope
        slope.windowidx = dsearchn(timebase_adj, makeCol(slope.window));
        slope.windowidx = [slope.windowidx(1) slope.windowidx(end)];

        movies = VSDmov.data(:,:,idxrange,:) ;

        %----------------------------------------------------------------
        % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
        %----------------------------------------------------------------
        rej = 'reject' ;
        if reject_on > 1
            rej = [rej num2str(reject_on)];
        end
        
        rejectidx = [];
        
        if setting.manual_reject
            try
                rejectidx = [rejectidx  makeRow(VSDI.(rej).manual)];
            catch
                rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
                disp(['rejec.manual was used for fish' num2str(VSDI.ref) 'because there is no reject' num2str(reject_on) '.manual'])
            end
        end
        
        if setting.GSabsthres_reject
            rejectidx = [rejectidx  makeRow(VSDI.(rej).GSabs025)];
            
        end
        
        if setting.GSmethod_reject
            rejectidx = [rejectidx makeRow(VSDI.(rej).GSdeviat2sd)];
        end
        
        if setting.force_include
            rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
            
        end
        
        rejectidx = sort(unique(rejectidx));
        
        % -------------------------------------------
        % LOOP THROUGH CONDITIONS (TO COMPUTE ONLY THE CONDITIONS OF
        % INTEREST)
        % -------------------------------------------
        
        for condition =  makeRow(cond_codes)
            % -------------------------------------------
            % SELECT CASES  AND AVERAGE MOVIE
            % -------------------------------------------
            sel_trials = find(VSDI.condition(:,1)==condition);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            % -------------------------------------------
            % CALCULATE MEASURES for each pixel and trial (and
            % store in 'maps')
            % -------------------------------------------
            for triali = makeRow(sel_trials)
                for xi = 1:size(movies,1)
                    for yi = 1:size(movies,2)
                        
                        pixelwave = squeeze(movies(xi, yi,:,triali));
                        pixelwave(end) = NaN;
                        temp = devo_peak2peak(pixelwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
                        
                        maps.peak(xi,yi,triali) = temp.peakminusbasel;
                        maps.wmean(xi,yi,triali) = temp.wmean;
                        
                        %slope
                        idx0 = slope.windowidx(1);
                        idxend = slope.windowidx(end);
                        waveW = pixelwave(idx0:idxend);
                        slopemean = mean(diff(waveW));
                        maps.slope(xi,yi,triali) = slopemean;
                        clear pixelwave waveW slopemean
                        
                    end % for yi
                end % for xi
                
                
            end % for triali
            blob()
        end %for condition
        
%         %% SPATIAL Z-SCORE AMONG CONDITIONS (for each measure)
%         dim = size(maps.peak);
%         temp_peak = reshape(maps.peak, [dim(1)*dim(2) dim(3)]);
%         temp_wmean = reshape(maps.wmean, [dim(1)*dim(2) dim(3)]);
%         temp_slope = reshape(maps.slope, [dim(1)*dim(2) dim(3)]);
%         
%         Zpeak = zscore(maps.peak, 0, 'all');
%         Zwmean = zscore(maps.wmean, 0, 'all');
%         Zslope = zscore(maps.slope, 0, 'all');
%         
%         mapsZ.peak =reshape(Zpeak, [dim(1) dim(2) dim(3)]);
%         
%         mapsZ.wmean =reshape(Zwmean, [dim(1) dim(2) dim(3)]);
%         
%         mapsZ.slope = reshape(Zslope, [dim(1) dim(2) dim(3)]);
  
%% PLOT 
maxval.peak = max(maps.peak(:))*0.9;
maxval.wmean = max(maps.wmean(:))*0.9; 
     
        for condition =  makeRow(cond_codes)
            % -------------------------------------------
            % SELECT CASES  AND AVERAGE MOVIE
            % -------------------------------------------
            sel_trials = find(VSDI.condition(:,1)==condition);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            
            % -------------------------------------------
            % CALCULATE MEASURES for each pixel and condition (and
            % store in 'maps')
            % -------------------------------------------
            figure
            for triali = makeRow(sel_trials)
                subplot(1,2,1);
                imagesc( maps.peak(:,:,triali) )
                axis image
                set(gca, 'clim', [0 maxval.peak])
                title('peak')
                colorbar
                
                subplot(1,2,2);
                imagesc( maps.wmean(:,:,triali) )
                axis image
                
                set(gca, 'clim', [0 maxval.wmean])
                title('wmean')
                colorbar
                
                sgtitle([num2str(VSDI.ref),'-cond', num2str(condition),'-trial idx',  num2str(triali)]);
                
                pause
            end % for triali
            
        end %for condition

     
            
%             % ------------------------------------------------------------------
%             ... PLOT MAPS OF SPATIAL Z-SCORE OF AVERAGE MEASURES CONDITION-WISE
%             % ------------------------------------------------------------------
%             
%             % colormap like jet but changing initial values 
%             ccmap = jet;
%             % remove darker colors: 
%             darkblue = ccmap(4,:);
%             ccmap = removerows(ccmap,'ind',[1:8, 1:3]);
% 
%             % change first color and interpolate
%             ccmap(1,:) = darkblue;
%             flag = 3;
%             R = linspace(ccmap(1,1), ccmap(flag,1), flag);
%             G = linspace(ccmap(1,2), ccmap(flag,2), flag);
%             B = linspace(ccmap(1,3), ccmap(flag,3), flag);
%             ccmap(1:flag,:) = [R; G; B]'; 
%             
%             
%             cclim = [0 4];
%             %         ccmap(1,:) = [0 0 0];
%             
% 
%             % PLOT PEAK + WMEAN
%             % ------------------------------------------------------------------
%             figure
%             ci = 0;
%             for condition =  makeRow(cond_lohi)
%                 % peak in the first row
%                 ci = ci+1;
%                 ax(ci) = subplot(2,ncond,ci);
%                 
%                 im1 = mapsZ.peak(:,:,ci);
%                 %                 im1 = interp2(im1, 5, 'nearest');
%                 %                 im1(~VSDI.crop.preview) = 0;
%                 
%                 alphamask = ones(size(im1))*0.6;
%                 
%                 %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
%                 %                 alphamask =  interp2(alphamask, 5, 'nearest');
%                 
%                 back = VSDI.backgr(:,:,VSDI.nonanidx(1));
%                 back = interp2(back,5);
%                 %                 imagesc(im1)
%                 %                 axis image
%                 %                 set(gca, 'clim',cclim)
%                 %                 colorbar; colormap(ccmap)
%                 %
%                 plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
%                 ax(ci).Visible = 'off';
%                 % STOPS THE CODE FOR CHECKING ROI CENTERS
%                 %             if ci ==2
%                 %                 return
%                 %             end
%                 
%                 condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
%                 tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
%                 title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
%                 
%                 % wmean in the second row
%                 ax(ci+ncond) = subplot(2,ncond,ci+ncond);
%                 
%                 im2 = mapsZ.wmean(:,:,ci);
%                 %                 im2(~VSDI.crop.preview) = 0;
%                 
%                 alphamask = ones(size(im2))*0.6;
%                 %                 alphamask= im2  >0;
%                 
%                 back = VSDI.backgr(:,:,VSDI.nonanidx(1));
%                 back = interp2(back,5);
%                 plot_framesoverlaid2(im2, back, alphamask, 0, ax(ci+ncond), cclim, 1 , 0, ccmap)
%                 ax(ci+ncond).Visible = 'off';
%                 
%             end
%             
%             sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Zmap. up: peak; down: onsetA. Cond:' num2str(cond_lohi)])
%             localname = ['MAPS_Zscore'] ;
%             
%             if savemaps
%                 
%                 %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
%                 set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
%                 print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
%                 
%                 close
%             end
            
        %
        
    
   blob()
    
% Last Update: 15/02/22