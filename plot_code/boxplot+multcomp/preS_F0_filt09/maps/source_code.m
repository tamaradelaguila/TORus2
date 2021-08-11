% Las medidas se aplican al average de cada condición

%% VISUALIZE MEASURES PIXEL-WISE.
clear

%----------------------------------------------------------------
% SELECT MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------

% measure = {'slopemean'};
% measure = {'slopemax'};
% measure ={'onsetnoise_ms'};
% measure = {'peakminusbasel'};
% measure = {'peakminusbasel' 'slopemean' 'slopemax' 'onsetnoise_ms' };
measure = {'peakminusbasel' 'slopemax' };

saveon = 0;

for nfish =  2%@ SET
    % user_settings
    path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/VSDI_ourToolbox/';
    path.TORus = '/home/tamara/Documents/MATLAB/VSDI/TORus';
    path.data = fullfile(path.TORus, 'data');
    path.grouplist = path.TORus;
    path.list =fullfile(path.TORus, 'data','BVlists');
    addpath(genpath(path.rootpath));
    addpath(path.TORus);
    % end of user_settings
    
    [VSDI] = TORus('load',nfish);
    
    temp = TORus('loadmovie',nfish,'_09filt3');
    movies = temp.data(:,:,1:end-1,:);
    
    % windows in which analyse the measures
    
    % windows = [0 100; 0 200; 100 200; 100 300] ; %@SET
    
    cond_codes = [200:203];
    
    %     cond_codes = unique(VSDI.condition(:,1));
    %     cond_codes=  cond_codes(~isnan(cond_codes));
    %     cond_codes = setdiff(cond_codes, 0);
    %     cond_codes =[2000 2001 2002 2003 3000 3001 3002 4000 4001 4002 ]; %for nfish = 6 (#210412)
    
    for reject_on =[0]  %@ SET
        
        setting.manual_reject = 0; %@ SET
        setting.GSmethod_reject = 1;  %@ SET
        setting.GSabsthres_reject = 1; %@ SET
        setting.force_include = 0; %@ SET
        
        % out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_measurements' ; %@ SET
        out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/boxplot/boxplot+multcomp/preS_F0_filt09/maps/intensos';
        %% SELECT CASES
        
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
        
        
        %% 1. PEAK-TO-PEAK MEASUREMENTS OF AVERAGE MOVIE
        
        % all windows for the function in 'ms'
        window.min = [-100 100];
        window.max = [0 600];
        window.movsum = 50;
        window.baseline = [-300 0];
        window.slope = 50;
        
        noise.SDfactor=4;
        noise.fr_abovenoise=30;
        method = 'movsum';
        
        %--------------------------------------
        %1. APPLY FUNCTION TO THE AVERAGE-MOVIE FOR EACH CONDITION AND PLOT
        %--------------------------------------
        j = 1;
        
        for condi = makeRow(cond_codes)
            sel_trials  = find(VSDI.condition(:,1)==condi);
            if reject_on  %@ SET
                sel_trials = setdiff(sel_trials, rejectidx);
            end
            
            avemovie = mean(movies(:,:,:,sel_trials),4);
            
            for rowi = 1:size(avemovie,1)
                for coli = 1:size(avemovie,2)
                    wave = squeeze(avemovie(rowi, coli, :));
                    output = devo_peak2peak(wave, VSDI.timebase, window,noise, method, 0);
                    
                    frames.peak2peak(rowi,coli,j) = output.p2p_value;
                    frames.peakminusbasel(rowi,coli,j) = output.peakminusbasel;
                    frames.peaklat(rowi,coli,j) = output.peaklat_ms;
                    frames.p2plat(rowi,coli,j) = output.p2plat_ms;
                    frames.onset30_latency_ms(rowi,coli,j) = output.onset30_latency_ms;
                    frames.onsetnoise_ms(rowi,coli,j) = output.onsetnoise_ms;
                    frames.noisethresh(rowi,coli,j) = output.noisethresh;
                    frames.slopemax(rowi,coli,j) = output.slopemax;
                    
                    peakidx.tmin(rowi,coli,j) = output.peakidx(1); %it'll be used for the calculation o
                    peakidx.tmax(rowi,coli,j) = output.peakidx(2); %it'll be used for the calculation o
                    
                    clear output
                end %coli
            end %rowi
            
            j = j+1;
            display(condi)
            clear sel_trials
            
        end %condi
        %         blob()
        
        %2. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        maxval.peak2peak= max(abs(frames.peak2peak(:)));
        maxval.peakminusbasel= max(abs(frames.peakminusbasel(:)));
        maxval.peaklat= max(abs(frames.peaklat(:)));
        maxval.p2plat= max(abs(frames.p2plat(:)));
        maxval.onset30_latency_ms= max(abs(frames.onset30_latency_ms(:)));
        maxval.onsetnoise_ms= max(abs(frames.onsetnoise_ms(:)));
        maxval.noisethresh= max(abs(frames.noisethresh(:)));
        maxval.slopemax= max(abs(frames.slopemax(:)));
        
        c_lim.peak2peak = [0 maxval.peak2peak(~isoutlier(maxval.peak2peak))];
        c_lim.peakminusbasel = [0 maxval.peakminusbasel(~isoutlier(maxval.peakminusbasel))];
        c_lim.peaklat = [0 maxval.peaklat(~isoutlier(maxval.peaklat))];
        c_lim.p2plat = [0 maxval.p2plat(~isoutlier(maxval.p2plat))];
        c_lim.onset30_latency_ms = [0 maxval.onset30_latency_ms];
        
        c_lim.onsetnoise_ms = [0 150];
        
        c_lim.noisethresh = [0 maxval.noisethresh(~isoutlier(maxval.noisethresh))];
        c_lim.slopemax = [0 maxval.slopemax(~isoutlier(maxval.slopemax))]; %smaller color limits
        
        BVmap = colormap_loadBV();
        
        
        for resulti = 1:length(measure)
            
            result= measure{resulti};
            %----------------------------------------------------------------
            % CONFIGURATION OF PARAMETERS THAT ARE SPECIFIC FOR EACH MEASURE
            %----------------------------------------------------------------
            measureframe= [];
            localclim = [];
            
            switch result
                
                % ------------------------------------------------------------------------
                case 'peakminusbasel'
                    
                    measureframe = frames.peakminusbasel;
                    
                    localmap = jet;
                    localclim = c_lim.peakminusbasel;
                    
                    
                    localtitle = [num2str(VSDI.ref), 'peak-b for each cond'];
                    localname = [num2str(VSDI.ref), '-peak-b(average mov)settings1'];
                    
                    
                    % ------------------------------------------------------------------------
                case 'slopemax'
                    
                    measureframe = frames.slopemax;
                    
                    localmap = jet;
                    localclim = c_lim.slopemax;
                    
                    localtitle = [num2str(VSDI.ref), 'slopemax for each cond'];
                    localname = [num2str(VSDI.ref), '-slopemax(average mov)settings1'];
                    
                    
                    % ------------------------------------------------------------------------
                case 'onsetnoise_ms'
                    
                    measureframe = frames.onsetnoise_ms;
                    
                    localmap = flipud(jet);
                    localclim = c_lim.onsetnoise_ms;
                    
                    localtitle = [num2str(VSDI.ref), 'onsetnoise_m_s for each cond'];
                    localname = [num2str(VSDI.ref), '-onsetnoise_ms(average mov)settings1'];
                    
                    % ------------------------------------------------------------------------
                    
                case 'noisethresh'
                    
                    measureframe = frames.noisethresh;
                    
                    localmap = jet;
                    localclim = c_lim.noisethresh;
                    
                    
                    localtitle = [num2str(VSDI.ref), 'peak-b for each cond'];
                    localname = [num2str(VSDI.ref), '-peak-b(average mov)settings1'];
                    
                    % ------------------------------------------------------------------------
                case 'slopemean' % in this case, the measure has to be computed
                    
                    % 1. CALCULATE MEASURE 'measureframe'
                    j = 1;
                    
                    for condi = makeRow(cond_codes)
                        sel_trials  = find(VSDI.condition(:,1)==condi);
                        if reject_on  %@ SET
                            sel_trials = setdiff(sel_trials, rejectidx);
                        end
                        
                        avemovie = mean(movies(:,:,:,sel_trials),4);
                        
                        for rowi = 1:size(avemovie,1)
                            for coli = 1:size(avemovie,2)
                                wave = squeeze(avemovie(rowi, coli, :));
                                output = devo_peak2peak(wave, VSDI.timebase, window, noise, method, 0);
                                
                                idx0= dsearchn(VSDI.timebase, 0);%get 0 index
                                waveW = wave(idx0:output.peakidx(2));
                                slopemean = mean(diff(waveW));
                                
                                measureframe(rowi, coli, j) = slopemean;
                            end %coli
                        end %rowi
                        j = j+1;
                        clear sel_trials
                    end % condi
                    
                    %GET MAX AND PLOT THEM WITH THE SAME LIMITS
                    localmax = max(abs(measureframe(:)));
                    
                    
                    % 2. CONFIGURE THE OTHER PARAMETERS
                    localmap = jet;
                    localclim = [0 localmax];
                    
                    localtitle = [num2str(VSDI.ref), 'slopemean for each cond'];
                    localname = [num2str(VSDI.ref), '-slopemean(average mov)settings1'];
                    
            end % result case selection
            
            % %--------------------------------------
            % %  PLOT MEASURES
            % %--------------------------------------
            
            figure
            
            for ploti = 1:length(cond_codes)
                
                ax(ploti) = subplot(3,4,ploti);
                imagesc(measureframe(:,:,ploti))
                set (ax(ploti), 'clim', localclim)
                colormap(localmap)
                condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
                tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
                title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
                
            end
            
            % plot colorbar
            ax(11) = subplot(3,4,11);
            colorbar
            set (ax(11), 'clim', localclim)
            colormap(localmap)
            
            %plot brain
            ax(12) = subplot(3,4,12);
            imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
            colormap(ax(12), bone)
            
            if  reject_on
                sgtitle([localtitle '(cl)'])
            else
                sgtitle(localtitle)
                
            end
            
            % save image and close @SET
            if  reject_on
                out.name =[localname, 'block' , num2str(cond_codes(1)), '(clean)'];
            else
                out.name = [localname, 'block' , num2str(cond_codes(1))];            
            end
            
%             saveas(gcf,fullfile(out.folder,[out.name '.jpg']),'jpg')
%             close
        end %result loop
        
        
    end %reject_one
    
    clearvars -except measure saveon output
    
end %nfish
blob() ; pause(0.1); blob()