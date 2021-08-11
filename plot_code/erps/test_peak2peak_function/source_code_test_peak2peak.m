
clear
settings_test1

nfish =10;%@ SET
[VSDI] = TORus('load',nfish);


VSDroiTS = TORus('loadwave',nfish);

% windows in which analyse the measures
setting.reject_on =1;  %@ SET

    setting.manual_reject = 0; %@ SET
    setting.GSmethod_reject = 1;  %@ SET
    setting.GSabsthres_reject = 1; %@ SET
    setting.force_include = 1; %@ SET


%% SELECT CASES
%     
%     cond_codes = unique(VSDI.condition(:,1));
%     cond_codes=  cond_codes(~isnan(cond_codes));
%     cond_codes = setdiff(cond_codes, 0);
    cond_codes =[2001 2002 2003];
    
    
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
    
    
    
    %% SETTINGS FOR PEAK-FINDER
    % DEBUG mi = 1; ci = 2002; triali = 69
    
    %% LOOP THROUGH CONDITIONS
    for mi = 1:2
        method = method_kind{mi}; 
        for ci = cond_codes
            trials = find(VSDI.condition(:,1) == ci);
            trials = setdiff(trials, rejectidx);
            
            for triali = makeRow(trials)
                wave = squeeze(VSDroiTS.filt306.data(:,nroi,triali));
                [output] = devo_peak2peak(wave, VSDI.timebase, wind, method, 1);

                name = strcat(num2str(VSDI.ref),method,'-cond',num2str(ci) , '-', VSDI.list(triali).Name(end-7:end-2),'.jpg');
                title(name)
                saveas(gcf, fullfile(pathsave, name), 'jpg')
                close
            end
        end
        
    end %method
    
        %% ERPs from each condition
        cond_codes = sort(unique(VSDI.condition(:,1)));
        cond_codes = cond_codes(~isnan(cond_codes));
        
    for mi = 1:2
        method = method_kind{mi}; 
        for ci = makeRow(cond_codes)
            trials = find(VSDI.condition(:,1) == ci);
            trials = setdiff(trials, rejectidx);
            
                wave = squeeze(mean(VSDroiTS.filt306.data(:,nroi,trials),3));
                [output] = devo_peak2peak(wave, VSDI.timebase, wind, method, 1);

                name = strcat(num2str(VSDI.ref),method,'-AVERAGE cond',num2str(ci) ,'.jpg');
                title(name)
                saveas(gcf, fullfile(pathsave, name), 'jpg')
                close
        end
        
    end %method