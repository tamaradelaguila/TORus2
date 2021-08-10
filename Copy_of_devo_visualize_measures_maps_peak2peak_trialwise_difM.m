%% El código es el mismo que devo_visualize_measures_maps_peak2peak
... excepto en la película a la que se le aplican las medidas: aquí se aplican a cada trial (al que se le ha restado la media de las películas control)
... todo el ploteado es el mismo

% get AVE-0


%%% DEBUG!! ACCORDING TO THE 'trialwise' SCRIPT! SEE LINE 82 (condition
%%% loop should be apply after, when calculating the means!!! otherwise all
%%% the trials are calculated for each conditions and then, when averaged,
%%% lower downs the measure


%% VISUALIZE MEASURES PIXEL-WISE.
clear
for nfish =  [12]%@ SET
    
    % clear
    user_settings
    
    [VSDI] = TORus('load',nfish);
    
    
    temp = TORus('loadmovie',nfish,'_06filt3');
    movies = temp.data(:,:,1:end-1,:);
    
            %% GET CONDITION CODES
%             cond_codes = unique(VSDI.condition(:,1));
%             cond_codes=  cond_codes(~isnan(cond_codes));
%             cond_codes = setdiff(cond_codes, 0); %delete code=0 if present
        cond_codes =[400:404 ]; %for nfish = 6 (#210412)

    for reject_on = [1]  %@ SET
        
        setting.manual_reject = 1; %@ SET
        setting.GSmethod_reject = 1;  %@ SET
        setting.GSabsthres_reject = 1; %@ SET
        setting.force_include = 1; %@ SET
        
        
        out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_measurements' ; %@ SET
        
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
        
        %% SELECT CASES
        sel_trials= [];
        
        for condi = cond_codes %to make sure that only the conditions of interest are computed)
            condtrials = makeCol(find(VSDI.condition(:,1)==condi));
            sel_trials  = [sel_trials; condtrials];
        end
        sel_trials= sort(sel_trials);
        
        if reject_on  %@ SET
            sel_trials = setdiff(sel_trials, rejectidx);
        end

        %% 1. PEAK-2-PEAK MEASUREMENTS OF EACH TRIAL MOVIE
        
        
        window.min = [-100 100];
        window.max = [0 600];
        window.movsum = 50;
        window.basel = [-100 0];
        
        method = 'movsum';
        
        %--------------------------------------
        %1. APPLY FUNCTION TO EACH MOVIE FROM EACH CONDITION AND PLOT
        %--------------------------------------
        tic
        
        for condi = makeRow(cond_codes) % we loop through conditions because we need to substract the control block, and also we ensure only needed trials are computed
            control_code = force0ending(condi); % control code
            
            trials_cond  = intersect(sel_trials, find(VSDI.condition(:,1)==condi));
            if reject_on  %@ SET
                trials_cond = setdiff(trials_cond, rejectidx);
            end
            
            sel_controls  = find(VSDI.condition(:,1)==control_code);
            if reject_on  %@ SET
                sel_controls = setdiff(sel_controls, rejectidx);
            end
            
            avecontrol = mean(movies(:,:,:,sel_controls),4);
            avecontrol = squeeze(avecontrol);
            
            for triali = makeRow(trials_cond)
                movtrial = squeeze(movies(:,:,:,triali));
                movie = movtrial - avecontrol;
                
                for rowi = 1:size(movie,1)
                    for coli = 1:size(movie,2)
                        wave = squeeze(movie(rowi, coli, :));
                        output = devo_peak2peak(wave, VSDI.timebase, window,[], method, 0);
                        
                        peak2peak(rowi,coli,triali) = output.p2p_value;
                        peakminusbasel(rowi,coli,triali) = output.peakminusbasel;
                        peaklat(rowi,coli,triali) = output.peaklat_ms;
                        p2plat(rowi,coli,triali) = output.p2plat_ms;
                        onset30_latency_ms(rowi,coli,triali) = output.onset30_latency_ms;
                        onsetnoise_ms(rowi,coli,triali) = output.onsetnoise_ms;
                        noisethresh(rowi,coli,triali) = output.noisethresh;
                                                
                        clear output wave
                        
                    end %coli
                end %rowi
                
                clear trials_cond
                
            end %triali
            clear sel_controls
            display(condi)
        end %condi
        
        t2 = toc
        %         blob()
        
        
        % AVERAGE THE MEASURES ACROSS TRIALS
        j = 1;
        for condi = makeRow(cond_codes)
            trials_cond  = intersect(sel_trials, find(VSDI.condition(:,1)==condi));
            
            frames.peak2peak(:,:,j) = squeeze(mean(peak2peak(:,:,trials_cond),3));
            frames.peakminusbasel(:,:,j) = squeeze(mean(peakminusbasel(:,:,trials_cond),3));
            frames.peaklat(:,:,j) = squeeze(mean(peaklat(:,:,trials_cond),3));
            frames.p2plat(:,:,j) = squeeze(mean(p2plat(:,:,trials_cond),3));
            frames.onset30_latency_ms(:,:,j) = squeeze(mean(onset30_latency_ms(:,:,trials_cond),3));
            frames.onsetnoise_ms(:,:,j) = squeeze(mean(onsetnoise_ms(:,:,trials_cond),3));
            frames.noisethresh(:,:,j) = squeeze(mean(noisethresh(:,:,trials_cond),3));
            
            j= j+1;
        end %condi
        
        %2. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        maxval.peak2peak= max(abs(frames.peak2peak(:)));
        maxval.peakminusbasel= max(abs(frames.peakminusbasel(:)));
        maxval.peaklat= max(abs(frames.peaklat(:)));
        maxval.p2plat= max(abs(frames.p2plat(:)));
        maxval.onset30_latency_ms= max(abs(frames.onset30_latency_ms(:)));
        maxval.onsetnoise_ms= max(abs(frames.onsetnoise_ms(:)));
        maxval.noisethresh= max(abs(frames.noisethresh(:)));

        
        c_lim.peak2peak = [0 maxval.peak2peak];
        c_lim.peakminusbasel = [0 maxval.peakminusbasel];
        c_lim.peaklat = [0 maxval.peaklat];
        c_lim.p2plat = [0 maxval.p2plat];
        c_lim.onset30_latency_ms = [0 maxval.onset30_latency_ms];
        c_lim.onsetnoise_ms = [0 maxval.onsetnoise_ms];
        c_lim.noisethresh = [0 maxval.noisethresh];
        
        BVmap = colormap_loadBV();
        
        
        %--------------------------------------
        %  PLOT peak2peak
        %--------------------------------------
        
        figure
        localmap = jet; 
        for ploti = 1:length(cond_codes)
            
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.peak2peak(:,:,ploti))
            set (ax(ploti), 'clim', c_lim.peak2peak)
            colormap(localmap)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.peak2peak)
        colormap(localmap)
        
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
        localmap= jet;
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.peakminusbasel(:,:,ploti))
            set (ax(ploti), 'clim', c_lim.peakminusbasel)
            colormap(localmap)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.peakminusbasel)
        colormap(localmap)
        
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

% %--------------------------------------
% %  PLOT onsetnoise_ms
% %--------------------------------------
  

        figure
        localmap = hot;
        for ploti = 1:length(cond_codes)
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.onsetnoise_ms(:,:,ploti))
            set (ax(ploti), 'clim', c_lim.peaklat)
            colormap(localmap)
                condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
                tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.onsetnoise_ms)
        colormap(localmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)    
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), 'onsetnoise_m_s (trial-0 mov) for each cond. (cl)'])
        else
            sgtitle([num2str(VSDI.ref), 'onsetnoise_m_s (trial-0 mov) for each cond.'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-onsetnoise_ms(trial-0 mov).settings1(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-onsetnoise_ms(trial-0 mov).settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close       
        
%         %--------------------------------------
%         %  PLOT peaklat
%         %--------------------------------------
%         
%         
%         figure
%         localmap= hot;
%         for ploti = 1:length(cond_codes)
%             ax(ploti) = subplot(3,4,ploti);
%             imagesc(frames.peaklat(:,:,ploti))
%             set (ax(ploti), 'clim', c_lim.peaklat)
%             colormap(localmap)
%             condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
%             tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
%             title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
%             
%         end
%         
%         % plot colorbar
%         ax(11) = subplot(3,4,11);
%         colorbar
%         set (ax(11), 'clim', c_lim.peaklat)
%         colormap(localmap)
%         
%         %plot brain
%         ax(12) = subplot(3,4,12);
%         imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
%         colormap(ax(12), bone)
%         
%         if  reject_on
%             sgtitle([num2str(VSDI.ref), 'peaklat (trial-0 mov) for each cond. (cl)'])
%         else
%             sgtitle([num2str(VSDI.ref), 'peaklat (trial-0 mov) for each cond.'])
%             
%         end
%         
%         % save image and close @SET
%         if  reject_on
%             out.name =[num2str(VSDI.ref), '-peaklat(trial-0 mov).settings1(clean)'];
%         else
%             out.name =[num2str(VSDI.ref), '-peaklat(trial-0 mov).settings1'];
%             
%         end
%         
%         saveas(gcf,fullfile(out.folder,out.name),'jpg')
%         close
%         
%         %--------------------------------------
%         %  PLOT p2plat
%         %--------------------------------------
%         
%         figure
%         localmap= hot;
%         for ploti = 1:length(cond_codes)
%             ax(ploti) = subplot(3,4,ploti);
%             imagesc(frames.p2plat(:,:,ploti))
%             set (ax(ploti), 'clim', c_lim.p2plat)
%             colormap(localmap)
%             condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
%             tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
%             title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
%             
%         end
%         
%         % plot colorbar
%         ax(11) = subplot(3,4,11);
%         colorbar
%         set (ax(11), 'clim', c_lim.p2plat)
%         colormap(localmap)
%         
%         %plot brain
%         ax(12) = subplot(3,4,12);
%         imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
%         colormap(ax(12), bone)
%         
%         if  reject_on
%             sgtitle([num2str(VSDI.ref), 'p2plat (trial-0 mov) for each cond. (cl)'])
%         else
%             sgtitle([num2str(VSDI.ref), 'p2plat (trial-0 mov) for each cond'])
%             
%         end
%         
%         % save image and close @SET
%         if  reject_on
%             out.name =[num2str(VSDI.ref), '-p2plat(trial-0 mov).settings1(clean)'];
%         else
%             out.name =[num2str(VSDI.ref), '-p2plat(trial-0 mov).settings1'];
%             
%         end
%         
%         saveas(gcf,fullfile(out.folder,out.name),'jpg')
%         close
%         
%         
%         %--------------------------------------
%         %  PLOT onset30_latency_ms
%         %--------------------------------------
%         
%         figure
%         localmap = hot;
%         for ploti = 1:length(cond_codes)
%             ax(ploti) = subplot(3,4,ploti);
%             imagesc(frames.onset30_latency_ms(:,:,ploti))
%             set (ax(ploti), 'clim', c_lim.onset30_latency_ms)
%             colormap(localmap)
%             condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
%             tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
%             title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
%             
%         end
%         
%         % plot colorbar
%         ax(11) = subplot(3,4,11);
%         colorbar
%         set (ax(11), 'clim', c_lim.onset30_latency_ms)
%         colormap(localmap)
%         
%         %plot brain
%         ax(12) = subplot(3,4,12);
%         imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
%         colormap(ax(12), bone)
%         
%         if  reject_on
%             sgtitle([num2str(VSDI.ref), '30%onset lat (trial-0 mov) for each cond (cl)'])
%         else
%             sgtitle([num2str(VSDI.ref), '30%onset lat (trial-0 mov) for each cond'])
%             
%         end
%         
%         % save image and close @SET
%         if  reject_on
%             out.name =[num2str(VSDI.ref), '-onset30lat(trial-0 mov).settings1 (clean)'];
%         else
%             out.name =[num2str(VSDI.ref), '-onset30lat(trial-0 mov).settings1'];
%             
%         end
%         
%         saveas(gcf,fullfile(out.folder,out.name),'jpg')
%         close
        
        
        
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
                movtrial = mean(movies(:,:,:,triali),4);
                movtrial = squeeze(movtrial);
                movie = movtrial - avecontrol;
                
                for rowi = 1:size(movie,1)
                    for coli = 1:size(movie,2)
                        
                        wave = squeeze(movie(rowi, coli, :));
                        output = devo_peak2peak(wave, VSDI.timebase, window,[], method, 0);
                        
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
        c_lim.meanslope = [0 maxval.meanslope];
        
        
        BVmap = colormap_loadBV();
        
        
        %--------------------------------------
        %  PLOT mean-slope
        %--------------------------------------
        
        figure
        localmap = jet;
        for ploti = 1:length(cond_codes)
            
            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.meanslope(:,:,ploti))
            set (ax(ploti), 'clim', c_lim.meanslope)
            colormap(localmap)
            condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
            tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.meanslope)
        colormap(localmap)
        
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
    end %reject_one
    clear
    
end %nfish
blob() ; pause(0.1); blob()