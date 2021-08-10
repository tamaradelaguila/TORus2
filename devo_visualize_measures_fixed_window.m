%% VISUALIZE MEASURES PIXEL-WISE.

clear
user_settings

nfish =6;%@ SET
[VSDI] = TORus('load',nfish);


temp = TORus('loadmovie',nfish,'_06filt3');
movies = temp.data(:,:,1:end-1,:);

% windows in which analyse the measures

% windows = [0 100; 0 200; 100 200; 100 300] ; %@SET
windows = [60 160] ; %@SET for more than one window, use different rows

setting.reject_on =1;  %@ SET

    setting.manual_reject = 0; %@ SET
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET
    setting.force_include = 1; %@ SET

    
out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/boxplot_03' ; %@ SET

%% SELECT CASES
    
    cond_codes = unique(VSDI.condition(:,1));
    cond_codes=  cond_codes(~isnan(cond_codes));
    cond_codes = setdiff(cond_codes, 0);
    % cond_codes =[100 101 102 103 300 301 302 303];
    
    
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
    
    
    %% MEASURE: MEAN ACTIVITY IN TIMEWINDOW
    %1. WINDOW        
        %1. WINDOW

        wind_ms = windows(wi,:); %window from the loop
        
        wind_idx = find_closest_timeidx(wind_ms, VSDI.timebase);
        wind_actualms =  VSDI.timebase(wind_idx);
        wind_actualms = [wind_actualms(1) wind_actualms(end)]; % keep only extreme values (safety code in case the code is changed and a range is selected)
               
        
        %2. GET AVERAGE MOVIE FOR EACH CONDITION AND MEASURE AND SAVE INTO A MATRIX
%         ncond= length(cond_codes);
        j = 1;
        
        for condi = makeRow(cond_codes)
            sel_trials  = find(VSDI.condition(:,1)==condi);
            if setting.reject_on  %@ SET      
                sel_trials = setdiff(sel_trials, rejectidx);     
            end
            
            frame_cond = mean(movies(:,:,wind_idx,sel_trials),4); frame_cond= squeeze(mean(frame_cond,3));
            frames_meanact(:,:,j) = frame_cond;
            j = j+1;
            clear sel_trials frames_cond
        end
        
        %3. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        frames = frames_meanact;
        maxval= max(abs(frames(:)));
        c_lim = [-maxval maxval];
        BVmap = colormap_loadBV();
        
        figure
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames(:,:,ploti))
            set (ax(ploti), 'clim', c_lim)
            colormap(BVmap)
            title(['cond=',num2str(cond_codes(ploti))])
            
        end
        
        % plot colorbar
        ax(11) = subplot(4,4,11);
        colorbar
        set (ax(11), 'clim', c_lim)
        colormap(BVmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)    
        
        if  setting.reject_on
            sgtitle([num2str(VSDI.ref), '-mean act in window',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'-for each cond (cl)'])
        else
            sgtitle([num2str(VSDI.ref), '-mean act in window',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'-for each cond'])
            
        end
        
        % save image and close @SET
        if  setting.reject_on
            out.name =[num2str(VSDI.ref), '-mean Act from',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2)),'(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-mean Act from',num2str(wind_actualms(1)),'to',num2str(wind_actualms(2))];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close
        
        
        clear frames frames_meanact maxval c_lim ax
        
        
        %% MEASURE: MEAN SLOPE IN TIMEWINDOW
        % windows = [0 100; 0 200; 100 200] ; %@SET
        % windows = windows';
        %
        %1. WINDOW
        % wind_ms = [202 502]; %@ SET
        
        % wind_idx = find_closest_timeidx(wind_ms, VSDI.timebase);
        % wind_actualms =  VSDI.timebase(wind_idx);
        % wind_actualms = [wind_actualms(1) wind_actualms(end)]; % keep only extreme values
        %
        %2. GET THE FRAGMENT OF THE MOVIE
        idx_range = wind_idx(1):wind_idx(end);
        movieW = movies(:,:,idx_range,:);
        moviesSlope = diff(movieW, [], 3);
        meanslope = squeeze(mean(moviesSlope,3));%for each trial;
        
        % GET AVERAGED FRAME FOR EACH CONDITION
%         ncond= length(cond_codes);
        j = 1;
        
        for condi = makeRow(cond_codes)
            
            sel_trials  = find(VSDI.condition(:,1)==condi);
            
            if setting.reject_on  %@ SET
                sel_trials = setdiff(sel_trials, rejectidx);    
            end
            
            frame_cond = mean(meanslope(:,:,sel_trials),3);
            frames_meanslope(:,:,j) = frame_cond;
            j = j+1;
            
        end % condi
        
        clear movieW movieSlope meanslope
        
        %3. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        frames = frames_meanslope;
        maxval= max(abs(frames(:)));
        c_lim = [-maxval maxval];
        BVmap = colormap_loadBV();
        
        figure
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(4,4,ploti);
            imagesc(frames(:,:,ploti))
            set (ax(ploti), 'clim', c_lim)
            colormap(BVmap)
            title(['cond=',num2str(cond_codes(ploti))])          
        end
        
                % plot colorbar
        ax(11) = subplot(4,4,11);
        colorbar
        set (ax(11), 'clim', c_lim)
        colormap(BVmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)    

        
        if  setting.reject_on
            sgtitle([num2str(VSDI.ref), '-mean slope in window',num2str(wind_actualms(1)),'to',num2str(wind_actualms(end)),'-for each cond (cl)'])
        else
            sgtitle([num2str(VSDI.ref), '-mean slope in window',num2str(wind_actualms(1)),'to',num2str(wind_actualms(end)),'-for each cond'])
            
        end
        
        % save image and close @SET
        if  setting.reject_on
            out.name =[num2str(VSDI.ref), '-mean slope from',num2str(wind_actualms(1)),'to',num2str(wind_actualms(end)),'(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-mean slope from',num2str(wind_actualms(1)),'to',num2str(wind_actualms(end))];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg') %out.folder defined in the header
        close
        
        clear frames frames_meanslope frames maxval c_lim ax
        wi
    end % wi
    
