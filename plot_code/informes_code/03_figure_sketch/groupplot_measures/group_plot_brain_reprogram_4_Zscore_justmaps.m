% code identical to:
... group_plot_brain_reprogram_4_Zscore2_working.m
... but to avoid having to mix groups


% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% BLANK SUBSTRACTION OF THE WAVES

% average movie
% extract wave for each pixel and make the measure maps (for all
% conditions): for condition for xi for yi

% z-spatial all conditions for each fish
% extract roi from z-spatially maps >>> definite measures
close all
clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus';
user_settings
cd(W)

load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new2.mat') % ATTENTION: select cases manually (there are repetitions)

% ///////////////////////////////////////////////////////////
% SETTINGS

selroinames = {'dm4m_R2','dm2_R2'}; %2roi

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_18filt6';% '_17filt5' ; '_18filt6'
% ref_movie= '_12filt5' ;

activ_unit = 'diffF'; % @ MANUALLY SET (just for display purposes)
analysisref = 'panel2_group8'; % MANUALLY SET!!! extra info for the name. Set group according to the rows of 'groupplot' selected in;   for suji  =  [1 3 4 9] lateral_group7_RECOV
reject_on = 3;

onset_factor = 0.1; % of max value to set rising threshold (to onset latency)

plotmaps = 1;
savemaps = 1; % if plotmaps =1

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/Zscore2' ;%@ SET

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms Range of analysis

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 1200]; % where to find the max peak
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

% END USER_SETTINGS
% ///////////////////////////////////////////////////////////////

% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{4,1} = 'wmean (ms)';
params{4,2} = feedf.window.wmean;

params{5,1} = 'slope win (ms)';
params{5,2} = slope.window;

params{6,1} = ['onset lat ' num2str(onset_factor*100) '%(ms)'];


if strcmpi( analysisref, 'panel2_group8') 
    sel_subjects = [1:4];
end

    %----------------------------------------------------------------
    % @SET: REJECT SETTINGS
    %----------------------------------------------------------------
    
    % Subsettings:
    setting.manual_reject = 1; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET+
    setting.force_include = 0; %@ SET
    

    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    i = 2; % counter for long format rows: the first will be the labels
    si = 0; % counter for subjects list
    
    for suji  =  sel_subjects %1:size(groupplot,1) SELECT included fish+condition
        
        si = si+1;
        nfish = groupplot{suji,1};
        cond_list = groupplot{suji,3};
        
        VSDI = TORus('load', nfish);
        VSDmov = TORus('loadmovie',nfish,ref_movie);
        movies = VSDmov.data ;
        F0 = VSDmov.F0;
        
        %----------------------------------------------------------------
        % GET INDEXES OF TIMERANGE
        %----------------------------------------------------------------
        idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
        idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values
        
        timebase_adj = VSDI.timebase(idxrange);
        
        % get indexes for slope
        slope.windowidx = dsearchn(timebase_adj, makeCol(slope.window));
        slope.windowidx = [slope.windowidx(1) slope.windowidx(end)];

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
                disp(['reject.manual was used for fish' num2str(VSDI.ref) 'because there is no reject' num2str(reject_on) '.manual'])
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
        
        
%         params{8,1} = 'rejected';

%         params{8,1} = rejectidx';

        
        % -------------------------------------------
        % FIRST LOOP THROUGH CONDITIONS AND PIXELS TO GET PIXEL-MEASURES
        % -------------------------------------------
        
        ci = 0;
        for condition =  makeRow(cond_list)
            ci = ci+1;
            
            cond_blank = force0ending(condition);
            % -------------------------------------------
            % SELECT CASES  AND AVERAGE MOVIE
            % -------------------------------------------
            sel_trials = find(VSDI.condition(:,1)==condition);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            % --------------------------------------------------------------------------
            % GET AVERAGE MOVIE. Use timerange set
            % --------------------------------------------------------------------------
                            
                    movieave = mean(movies(:,:,idxrange,sel_trials),4);
                    movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    ...(that normally corresponds to the background)
                            
            % -------------------------------------------
            % CALCULATE MEASURES for each pixel and condition (and
            % store in 'maps')
            % -------------------------------------------
                    
                    for xi = 1:size(movieave,1)
                        for yi = 1:size(movieave,2)
                            pixelwave = movieave(xi, yi,:);
                            
                            temp = devo_peak2peak(pixelwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
                            
                            maps.peak(xi,yi,ci) = temp.peakminusbasel;
                            maps.wmean(xi,yi,ci) = temp.wmean;
%                             maps.peaklat(xi,yi,ci) = temp.peaklatency;
%                             % pixel-wise peak latency
                            
                            %slope
                                idx0 = slope.windowidx(1);
                                idxend = slope.windowidx(end);
                                waveW = pixelwave(idx0:idxend);
                                slopemean = mean(diff(waveW));
                                maps.slope(xi,yi,ci) = slopemean;
                                clear pixelwave waveW slopemean

                        end % for yi
                    end % for xi

        end % for condition


        %% SPATIAL Z-SCORE AMONG CONDITIONS (for each measure)
        
%         dim = size(maps.peak);
%         temp_peak = reshape(maps.peak, [dim(1)*dim(2) dim(3)]);
%         temp_wmean = reshape(maps.wmean, [dim(1)*dim(2) dim(3)]);
%         temp_slope = reshape(maps.slope, [dim(1)*dim(2) dim(3)]);
%         temp_onsetlat = reshape(maps.slope, [dim(1)*dim(2) dim(3)]);
        
        Zpeak = zscore(maps.peak, 0, 'all');
        Zwmean = zscore(maps.wmean, 0, 'all');
        Zslope = zscore(maps.slope, 0, 'all');
        
        mapsZ.peak = Zpeak;
        mapsZ.wmean = Zwmean;
        mapsZ.slope = Zslope;

%         mapsZ.peak =reshape(Zpeak, [dim(1) dim(2) dim(3)]);
%         mapsZ.wmean =reshape(Zwmean, [dim(1) dim(2) dim(3)]);
%         mapsZ.slope = reshape(Zslope, [dim(1) dim(2) dim(3)]);

  
        
        %% PLOT MAPS
        if plotmaps
            
            savename = [num2str(VSDI.ref) '_' ref_movie '_rej' num2str(reject_on) '_'];
            %
            ncond = length(cond_list);
            
            %             % ------------------------------------------------------------------
            %             % PLOT MAPS OF AVERAGE MEASURES CONDITION-WISE
            %             % ------------------------------------------------------------------
            %             figure
            %             ci = 0;
            %             for condition =  makeRow(cond_list)
            %                 % peak in the first row
            %                 ci = ci+1;
            %                 ax(ci) = subplot(2,ncond,ci);
            %                 imagesc(maps.peak(:,:,ci))
            %                 axis image
            %                 colorbar
            %                 %
            %                 set(gca, 'clim', [0 max(maps.peak(:))*0.8])
            %                 condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
            %                 tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            %                 title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
            %
            %                 % wmean in the second row
            %                 ax(ci+ncond) = subplot(2,ncond,ci+ncond);
            %                 imagesc(maps.wmean(:,:,ci))
            %                 axis image
            %                 colorbar
            %                 set(gca, 'clim', [0 max(maps.wmean(:))*0.8])
            %
            %
            %             end
            %
            %             sgtitle([num2str(VSDI.ref), '(',ref_movie,  ')', 'up: peak; down: onsetA'])
            %             localname = ['MAPS'] ;
            %
            %             if savemaps
            %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
            %                 close
            %             end
            
%             % ------------------------------------------------------------------
%             ... PLOT MAPS OF SPATIAL Z-SCORE OF AVERAGE MEASURES CONDITION-WISE
%             % ------------------------------------------------------------------
            
            % colormap like jet but changing initial values 
            ccmap = jet;
            % remove darker colors: 
            darkblue = ccmap(4,:);
            ccmap = removerows(ccmap,'ind',[1:8, 1:3]);

            % change first color and interpolate
            ccmap(1,:) = darkblue;
            flag = 3;
            R = linspace(ccmap(1,1), ccmap(flag,1), flag);
            G = linspace(ccmap(1,2), ccmap(flag,2), flag);
            B = linspace(ccmap(1,3), ccmap(flag,3), flag);
            ccmap(1:flag,:) = [R; G; B]'; 
            
            cclim = [0 4];
            %         ccmap(1,:) = [0 0 0];
            

            % PLOT PEAK + WMEAN
            % ------------------------------------------------------------------
            figure
            ci = 0;
            for condition =  makeRow(cond_list)
                % peak in the first row
                ci = ci+1;
                ax(ci) = subplot(2,ncond,ci);
                
                im1 = mapsZ.peak(:,:,ci);
                %                 im1 = interp2(im1, 5, 'nearest');
                %                 im1(~VSDI.crop.preview) = 0;
                
                alphamask = ones(size(im1))*0.6;
                
                %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
                %                 alphamask =  interp2(alphamask, 5, 'nearest');
                
                back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                back = interp2(back,5);
                %                 imagesc(im1)
                %                 axis image
                %                 set(gca, 'clim',cclim)
                %                 colorbar; colormap(ccmap)
                %
                plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
                ax(ci).Visible = 'off';
                % STOPS THE CODE FOR CHECKING ROI CENTERS
                %             if ci ==2
                %                 return
                %             end
                
                condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
                tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
                title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
                
                % wmean in the second row
                ax(ci+ncond) = subplot(2,ncond,ci+ncond);
                
                im2 = mapsZ.wmean(:,:,ci);
                %                 im2(~VSDI.crop.preview) = 0;
                
                alphamask = ones(size(im2))*0.6;
                %                 alphamask= im2  >0;
                
                back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                back = interp2(back,5);
                plot_framesoverlaid2(im2, back, alphamask, 0, ax(ci+ncond), cclim, 1 , 0, ccmap)
                ax(ci+ncond).Visible = 'off';
                
            end
            
            sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Zmap. up: peak; down: onsetA. Cond:' num2str(cond_list)])
            localname = ['MAPS_Zscore2'] ;
            
            if savemaps
                
                %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
                
                close
            end
            

            % PLOT SLOPE
            % ------------------------------------------------------------------
            figure
            ci = 0;
            for condition =  makeRow(cond_list)
                % peak in the first row
                ci = ci+1;
                ax(ci) = subplot(1,ncond,ci);
                
                im1 = mapsZ.slope(:,:,ci);
                %                 im1 = interp2(im1, 5, 'nearest');
                %                 im1(~VSDI.crop.preview) = 0;
                
                alphamask = ones(size(im1))*0.6;
                
                %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
                %                 alphamask =  interp2(alphamask, 5, 'nearest');
                
                back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                back = interp2(back,5);
                %                 imagesc(im1)
                %                 axis image
                %                 set(gca, 'clim',cclim)
                %                 colorbar; colormap(ccmap)
                %
                plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
                ax(ci).Visible = 'off';
                % STOPS THE CODE FOR CHECKING ROI CENTERS
                %             if ci ==2
                %                 return
                %             end
                
            end
            
            sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Zmap. Slope. Cond:' num2str(cond_list)])
            localname = ['Zscore2_slope'] ;
            
            if savemaps
                
                %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
                
                close
            end
            
            
        end % if plotmaps
        
        
        %         % PLOT PEAK VALUES
        %         ax2 = subplot(1,3,2);
        %         temp = barplot.peak(suji,:,selroi)';
        %         legend(selroinames)
        %         title(num2str(suji))
        %
        %         localref = [num2str(VSDI.ref) '- cond=' num2str(condition) '-' ref_movie '-' refcase  '-reject' num2str(reject_on) ];
        %         sgtitle(localref)
        %
        %         savename = [ localref   '-' num2str(numel(selroi)) roikind 'ROI - PART2 WAVES'];
        % %         localname = 'tilewaves'];
        %         saveas(gcf, fullfile(savein, [savename localname '.jpg']), 'jpg')
        %         close

        clear maps mapsZ
        
        groupplot_print(suji,:) = groupplot(suji,1:end);
    end % for suji
    
    
    blob()

% Created 24/02/2022 from: 'group_plot_brain_reprogram_4_Zscore_working.m'
% Only extracts the maps, so the roi are not needed
