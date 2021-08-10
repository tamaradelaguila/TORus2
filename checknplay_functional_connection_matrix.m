nfish= 12;

VSDI = TORus('load', nfish);
VSDroiTS = TORus('loadwave',nfish) ;
waves = VSDroiTS.circ_filt309.data;

cond_codes = [400:404];
reject_on = 1;

setting.manual_reject = 0; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET
setting.force_include = 1; %@ SET



% SELECT EXCLUDED

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


%% CONCATENATING
% method = 'concatenate';
method = 'trialmean';

ms = [0;600];
msidx = dsearchn(VSDI.timebase, ms);

Rroi= [1 3 5 7 9 11 13];
Lroi= [2 4 6 8 10 12 14];

rois=  Rroi;
Rlabels = [VSDI.roi.labels(rois)];
cond_codes = [400 401 402 403 404];

j= 1;
figure
for condi = makeRow(cond_codes)
    trialidx=[];
    trialidx = setdiff(find(VSDI.condition(:,1) ==condi) , rejectidx);
    
    subplot(1,length(cond_codes),j)
    longwave = [];
    for tri = 1:length(trialidx)
        localwave = waves(msidx(1):msidx(2),rois,trialidx(tri));
        longwave = cat(1, longwave, localwave);
    end
    
    n_roi = size(longwave,2);
    
    switch method
        
        case 'concatenate'
            fcM= [];
            for roi1 = 1:n_roi
                for roi2 = 1:n_roi
                    % adjM(roi1,roi2, trialidx) = corr (sel_data(1:680,roi1,trialidx), sel_data(1:680,roi2,trialidx));
                    fcM (roi1,roi2) = corr(longwave(:,roi1), longwave(:,roi2));
                end % roi1
            end % roi2
            
        case 'trialmean'
                        fcM= [];

            for tri = 1:length(trialidx)
                for roi1 = 1:n_roi
                    for roi2 = 1:n_roi
                        wave1 =  waves(msidx(1):msidx(2),roi1,trialidx(tri));
                        wave2 =  waves(msidx(1):msidx(2),roi2,trialidx(tri));
                        % adjM(roi1,roi2, trialidx) = corr (sel_data(1:680,roi1,trialidx), sel_data(1:680,roi2,trialidx));
                        fcM (roi1,roi2,tri) = corr(wave1, wave2);
                        clear wave1 wave2
                    end % roi1
                end % roi2
                
            end %tri
            fcM = mean(fcM,3);
            
    end %switch
    
    imagesc(fcM); colorbar; set(gca, 'clim', [0, 1]   )
    %     plot(waves(:,3,trialidx(1))); hold; plot(longwave(1:340,3)); legend wave longwave
    title(['c',num2str(condi)])
    xticks(1:numel(rois)); xticklabels(Rlabels); xtickangle(90)
    yticks(1:numel(rois)); yticklabels(Rlabels);
    
    j=j+1;
end

sgtitle([num2str(VSDI.ref),'w:',num2str(ms'),'ms (' ,  method, ')' ])