%% El código es el mismo que devo_visualize_measures_maps_peak2peak
... excepto en la película a la que se le aplican las medidas: aquí se aplican a cada trial (al que se le ha restado la media de las películas control)
... todo el ploteado es el mismo

% get AVE-0
%% VISUALIZE MEASURES PIXEL-WISE.

 nfish =  10%@ SET
    
    % clear
    user_settings
    
    [VSDI] = TORus('load',nfish);
    
    
    temp = TORus('loadmovie',nfish,'_06filt3');
    movies = temp.data(:,:,1:end-1,:);
    
    reject_on = 1  %@ SET
        
        setting.manual_reject = 1; %@ SET
        setting.GSmethod_reject = 1;  %@ SET
        setting.GSabsthres_reject = 1; %@ SET
        setting.force_include = 1; %@ SET
        
        
        out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_measurements/test210509' ; %@ SET
        
        %% CONFIG REJECTION OPTIONS
        
        
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
        
        %% GET CONDITION CODES
            cond_codes = unique(VSDI.condition(:,1));
            cond_codes=  cond_codes(~isnan(cond_codes));
            cond_codes = setdiff(cond_codes, 0); %delete code=0 if present
%         cond_codes =[1001 1002 1003 2002 2003]; 
       
        %% 1. PEAK-2-PEAK MEASUREMENTS OF AVERAGE MOVIE
        

        window.min = [-100 100];
        window.max = [0 600];
        window.movsum = 50;
        window.basel = [-100 0];
        
        method = 'movsum';
        
        %--------------------------------------
        %1. APPLY FUNCTION TO THE AVERAGE-MOVIE FOR EACH CONDITION AND PLOT
        %--------------------------------------
        tic
        j = 1;
        
        for condi = makeRow(cond_codes)
            control_code = force0ending(condi); % control code
            
            sel_trials  = find(VSDI.condition(:,1)==condi);
            if reject_on  %@ SET
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            sel_controls  = find(VSDI.condition(:,1)==control_code);
            if reject_on  %@ SET
                sel_controls = setdiff(sel_controls, rejectidx);
            end
            
            avecontrol = mean(movies(:,:,:,sel_controls),4);
            avecontrol = squeeze(avecontrol);
            
            for triali = sel_trials
                movtrial = mean(movies(:,:,:,sel_trials),4);
                movtrial = squeeze(movtrial);
                movie = movtrial - avecontrol;
                
                for rowi = 1:size(movie,1)
                    for coli = 1:size(movie,2)
                        wave = squeeze(movie(rowi, coli, :));
                        output = devo_peak2peak(wave, VSDI.timebase, window, method, 0);
                        
                        frames.peak2peak(rowi,coli,triali,j) = output.p2p_value;
                        frames.peakminusbasel(rowi,coli,triali,j) = output.peakminusbasel;
                        frames.peaklat(rowi,coli,triali,j) = output.peaklat_ms;
                        frames.p2plat(rowi,coli,triali,j) = output.p2plat_ms;
                        frames.onset30_latency_ms(rowi,coli,triali,j) = output.onset30_latency_ms;
                        
                        peakidx.tmin(rowi,coli,triali,j) = output.peakidx(1); %it'll be used for the calculation o
                        peakidx.tmax(rowi,coli,triali,j) = output.peakidx(2); %it'll be used for the calculation o
                        
                        clear output
                        
                    end %coli
                end %rowi
                
                clear sel_trials
                
            end %triali
            clear sel_controls
            j = j+1;
            display(condi)
        end %condi
        
        t2 = toc
        %         blob()
        
        
        % AVERAGE ACROSS TRIALS
        frames.peak2peak = squeeze(mean(frames.peak2peak,3));
        frames.peakminusbasel = squeeze(mean(frames.peakminusbasel,3));
        frames.peaklat = squeeze(mean(frames.peaklat,3));
        frames.p2plat = squeeze(mean(frames.p2plat,3));
        frames.onset30_latency_ms = squeeze(mean(frames.onset30_latency_ms,3));
        
        peakidx.tmin = squeeze(mean(peakidx.tmin,3));
        peakidx.tmax = squeeze(mean(peakidx.tmax,3));
        
        
        %2. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        maxval.peak2peak= max(abs(frames.peak2peak(:)));
        maxval.peakminusbasel= max(abs(frames.peakminusbasel(:)));
        maxval.peaklat= max(abs(frames.peaklat(:)));
        maxval.p2plat= max(abs(frames.p2plat(:)));
        maxval.onset30_latency_ms= max(abs(frames.onset30_latency_ms(:)));
        
        c_lim.peak2peak = [-maxval.peak2peak maxval.peak2peak];
        c_lim.peakminusbasel = [-maxval.peakminusbasel maxval.peakminusbasel];
        c_lim.peaklat = [0 maxval.peaklat];
        c_lim.p2plat = [0 maxval.p2plat];
        c_lim.onset30_latency_ms = [0 maxval.onset30_latency_ms];
        
        BVmap = colormap_loadBV();
        

        %% PLOTS
        %--------------------------------------
        %  PLOT peak2peak
        %--------------------------------------
        
        figure
        climlocal = c_lim.peak2peak;
        for ploti = 1:length(cond_codes)
            
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.peak2peak(:,:,ploti))
            set (ax(ploti), 'clim', climlocal)
            colormap(BVmap)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', climlocal)
        colormap(BVmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), 'peak2peak (trial-0 mov) for each cond (cl)'])
        else
            sgtitle([num2str(VSDI.ref), 'peak2peak (trial-0 mov) for each cond'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-peak2peak(trial-0 mov).settings1(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-peak2peak(trial-0 mov).settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close
        
        
        %--------------------------------------
        %  PLOT peakminusbasel
        %--------------------------------------
        
        
        figure
        climlocal= c_lim.peakminusbasel * 0.3;
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.peakminusbasel(:,:,ploti))
            set (ax(ploti), 'clim', climlocal)
            colormap(BVmap)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', climlocal)
        colormap(BVmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), 'peak-basel (trial-0 mov) for each cond. (cl)'])
        else
            sgtitle([num2str(VSDI.ref), 'peak-basel (trial-0 mov) for each cond.'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-peakminusbasel(trial-0 mov).settings1(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-peakminusbasel(trial-0 mov).settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close
        
        %--------------------------------------
        %  PLOT peaklat
        %--------------------------------------
        
        % noise-threshold map
        plot_peaklat = frames.peaklat;
        noise_thresh = return_noisethresh(movie, VSDI.timebase, [-200 0], 2);
        outidx = frames.peakminusbasel < noise_thresh.sup & frames.peakminusbasel > noise_thresh.inf;
        plot_peakminusb = frames.peakminusbasel(outidx); 
        plot_peaklat(outidx) = maxval.peaklat;

        % plot
        figure
        climlocal = c_lim.peaklat;
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.peaklat(:,:,ploti))
            set (ax(ploti), 'clim', climlocal)
            colormap(flipud(jet))
            colorbar
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', climlocal)
        colormap(flipud(jet))
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), 'peaklat (trial-0 mov) for each cond. (cl)'])
        else
            sgtitle([num2str(VSDI.ref), 'peaklat (trial-0 mov) for each cond.'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-peaklat(trial-0 mov).settings1(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-peaklat(trial-0 mov).settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close
        
        %--------------------------------------
        %  PLOT p2plat
        %--------------------------------------
        
        figure
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.p2plat(:,:,ploti))
            set (ax(ploti), 'clim', c_lim.p2plat)
            colormap(jet)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.p2plat)
        colormap(jet)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), 'p2plat (trial-0 mov) for each cond. (cl)'])
        else
            sgtitle([num2str(VSDI.ref), 'p2plat (trial-0 mov) for each cond'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-p2plat(trial-0 mov).settings1(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-p2plat(trial-0 mov).settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close
        
        
        %--------------------------------------
        %  PLOT onset30_latency_ms
        %--------------------------------------
        
        figure
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.onset30_latency_ms(:,:,ploti))
            set (ax(ploti), 'clim', c_lim.onset30_latency_ms)
            colormap(jet)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.onset30_latency_ms)
        colormap(jet)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), '30%onset lat (trial-0 mov) for each cond (cl)'])
        else
            sgtitle([num2str(VSDI.ref), '30%onset lat (trial-0 mov) for each cond'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-onset30lat(trial-0 mov).settings1 (clean)'];
        else
            out.name =[num2str(VSDI.ref), '-onset30lat(trial-0 mov).settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close
        
        
        
        %% 2. PEAK-to-PEAK MEAN SLOPE OF AVERAGE MOVIE
        
        
        window.min = [-100 100];
        window.max = [0 600];
        window.movsum = 50;
        
        method = 'movsum';
        
        %--------------------------------------
        %1. APPLY FUNCTION TO THE AVERAGE-MOVIE FOR EACH CONDITION AND CALCULATE
        %MEAN SLOPE
        %--------------------------------------
        tic
        j = 1;
        
        for condi = makeRow(cond_codes)
            control_code = force0ending(condi); % control code
            
            sel_trials  = find(VSDI.condition(:,1)==condi);
            if reject_on  %@ SET
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            sel_controls  = find(VSDI.condition(:,1)==control_code);
            if reject_on  %@ SET
                sel_controls = setdiff(sel_controls, rejectidx);
            end
            
            avecontrol = mean(movies(:,:,:,sel_controls),4);
            avecontrol = squeeze(avecontrol);
            
            for triali = sel_trials
                movtrial = mean(movies(:,:,:,sel_trials),4);
                movtrial = squeeze(movtrial);
                movie = movtrial - avecontrol;
                
                for rowi = 1:size(movie,1)
                    for coli = 1:size(movie,2)
                        
                        wave = squeeze(movie(rowi, coli, :));
                        output = devo_peak2peak(wave, VSDI.timebase, window, method, 0);
                        
                        waveW = wave(output.peakidx(1):output.peakidx(2));
                        waveslope = diff(waveW);
                        meanslope = mean(waveslope);
                        
                        frames.meanslope(rowi, coli,triali, j) = meanslope;
                    end %COLI
                end %rowi
                clear sel_trials
                
            end %triali
            j = j+1
            clear sel_controls
            
        end %condi
        
        toc
        
        frames.meanslope = squeeze(mean(frames.meanslope,3));
        
        %2. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        maxval.meanslope= max(abs(frames.meanslope(:)));
        c_lim.meanslope = [-maxval.meanslope maxval.meanslope];
        
        
        BVmap = colormap_loadBV();
        
        
        %--------------------------------------
        %  PLOT mean-slope
        %--------------------------------------
        
        figure
        for ploti = 1:length(cond_codes)
            
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.meanslope(:,:,ploti))
            set (ax(ploti), 'clim', c_lim.meanslope)
            colormap(BVmap)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.meanslope)
        colormap(BVmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), 'mean slope in p2p (trial-0 mov) for each cond (cl)'])
        else
            sgtitle([num2str(VSDI.ref), 'mean slope in p2p (trial-0 mov) for each cond'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-p2p mean slope (trial-0 mov). settings1 (clean)'];
        else
            out.name =[num2str(VSDI.ref), '-p2p mean slope (trial-0 mov). settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close
        
        blob()
        
