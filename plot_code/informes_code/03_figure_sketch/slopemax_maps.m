%% El código es el mismo que devo_visualize_measures_maps_peak2peak, pero conservando sólo slopemax
... excepto en la película a la que se le aplican las medidas: aquí se aplican a la pelíucla average de cada condición a la que se le ha sustraído el average del cotrol 

... todo el ploteado es el mismo
    
% clear

 nfish = [11]%  [2 3 4 8 9 10]%@ SET

    % clear
    user_settings

    [VSDI] = TORus('load',nfish);


    temp = TORus('loadmovie',nfish,'_12filt5');
    movies = temp.data(:,:,1:end-1,:);


% GET CONDITION CODES   
%     cond_codes = unique(VSDI.condition(:,1));
%     cond_codes=  cond_codes(~isnan(cond_codes));
%     cond_codes = setdiff(cond_codes, 0); %delete code=0 if present
    cond_codes =[400:404]; %for nfish = 6 (#210412)

 reject_on =  [0]  %@ SET

    setting.manual_reject = 1; %@ SET
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET
    setting.force_include = 1; %@ SET

    
out.folder = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/test_measurements' ; %@ SET

% CONFIG REJECTION OPTIONS
    
    
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
    

    % 1. PEAK-2-PEAK MEASUREMENTS OF AVERAGE MOVIE
    
    window.min = [-100 100]; 
    window.max = [0 600];
    window.movsum = 50;
    window.basel = [-100 0];
    window.slope=50;
    window.wmean=[0 250];

    method = 'movsum';

%% SLOPEMAX FROM AVERAGE MOVIE

%--------------------------------------
%1. APPLY FUNCTION TO THE AVERAGE-MOVIE (minus Average control) FOR EACH CONDITION AND PLOT
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
            
            avetrials = mean(movies(:,:,:,sel_trials),4);
            avecontrol = mean(movies(:,:,:,sel_controls),4);
            avemovie = avetrials - avecontrol; 
            
            for rowi = 1:size(avemovie,1)
                for coli = 1:size(avemovie,2)
                wave = squeeze(avemovie(rowi, coli, :));
                output = devo_peak2peak(wave, VSDI.timebase, window, [], method, 0);
                
                frames.slopemax(rowi,coli,j) = output.slopemax;


                clear output
                end %coli
            end %rowi
          
            j = j+1;
            clear sel_trials 
            display(condi)
        end %condi
        
        t2 = toc
        blob()
        
        %2. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        maxval.slopemax= max(abs(frames.slopemax(:)));

        c_lim.slopemax = [-maxval.slopemax maxval.slopemax];

        BVmap = colormap_loadBV();
        
%--------------------------------------
%  PLOT slopemax
%--------------------------------------

        figure
        localmap = jet;
        for ploti = 1:length(cond_codes)

            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.slopemax(:,:,ploti))
%             set (ax(ploti), 'clim', c_lim.slopemax)
            colormap(localmap)
            colorbar
                condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
                tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.slopemax)
        colormap(localmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        colormap(ax(12), bone)    
        
        if  reject_on
            sgtitle([num2str(VSDI.ref), 'slopemax (avedif mov) for each cond (cl)'])
        else
            sgtitle([num2str(VSDI.ref), 'slopemax (avedif mov) for each cond'])
            
        end
        
        % save image and close @SET
        if  reject_on
            out.name =[num2str(VSDI.ref), '-slopemax(avedif mov).settings1(clean)'];
        else
            out.name =[num2str(VSDI.ref), '-slopemax(avedif mov).settings1'];
            
        end
        
        saveas(gcf,fullfile(out.folder,out.name),'jpg')
        close

blob() ; pause(0.1); 


%% SLOPEMAX TO EACH TRIAL AND LATER AVERAGE
   %--------------------------------------
        % APPLY FUNCTION TO EACH TRIAL AND PLOT
        %--------------------------------------
        
        
        for triali = makeRow(VSDI.nonanidx) % we loop through conditions because we need to substract the control block, and also we ensure only needed trials are computed
            condi = VSDI.condition(triali,1);
            control_code = force0ending(condi); % control code
            
                        
            movtrial = squeeze(movies(:,:,:,triali));
            
            for rowi = 1:size(movtrial,1)
                for coli = 1:size(movtrial,2)
                    wave = squeeze(movtrial(rowi, coli, :));
                    output = devo_peak2peak(wave, VSDI.timebase, window,[], method, 0, 0);
                    
                    slopemax_trials(rowi,coli,triali) = output.slopemax;
                    
                    clear output wave
                    
                end %coli
            end %rowi
            
            clear sel_controls
            display(triali)
            
        end %triali
        
        %         blob()
        
        %--------------------------------------
        % AVERAGE THE MEASURES CONDITION-WISE ACROSS TRIALS
        %--------------------------------------
        
        j = 1;
        for condi = makeRow(cond_codes)
            trials_cond  = intersect(sel_trials, find(VSDI.condition(:,1)==condi));
            avecond = squeeze(mean(slopemax_trials(:,:,trials_cond),3));
            
            cond0 = force0ending(condi);
            trials0 =  intersect(sel_trials, find(VSDI.condition(:,1)==cond0));
            avecontrol = squeeze(mean(slopemax_trials(:,:,trials0),3));
            
            frames.slopemax_trials(:,:,j) = avecond -avecontrol ;
            
            j= j+1
        end %condi
        
        %2. GET MAX-MIN AND PLOT THEM WITH THE SAME LIMITS
        maxval.slopemax_trials= max(abs(frames.slopemax_trials(:)));
        
        
        c_lim.noisethresh = [0 maxval.slopemax_trials];
        
        BVmap = colormap_loadBV();
        
        
        
                figure
        localmap = jet;
        for ploti = 1:length(cond_codes)

            ax(ploti) = subplot(3,4,ploti);
            imagesc(frames.slopemax_trials(:,:,ploti))
            axis image
%             set (ax(ploti), 'clim', c_lim.slopemax)
            colormap(localmap)
            colorbar
                condidx = find(VSDI.condition(:,1) ==cond_codes(ploti)); % get idx from condition
                tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            title(['c',num2str(cond_codes(ploti)), '(', num2str(tempmA),'mA)'])
            
        end
        
        % plot colorbar
        ax(11) = subplot(3,4,11);
        colorbar
        set (ax(11), 'clim', c_lim.slopemax_trials)
        colormap(localmap)
        
        %plot brain
        ax(12) = subplot(3,4,12);
        imagesc(VSDI.backgr(:,:,VSDI.nonanidx(1)));
        axis image
        colormap(ax(12), bone)    

        
   