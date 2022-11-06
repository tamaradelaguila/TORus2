clear
W = pwd;
cd '/home/tamara/Documents/MATLAB/VSDI/TORus'
user_settings
cd(W)

%----------------------------------------------------------------
% @SET: fish + conditions
%---------------------------------------------------------------
nfish = 22;
movie_ref = '_18filt6'; % input movie '_17filt5'

%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus ('load', nfish);
VSDmov = TORus('loadmovie',nfish,movie_ref);

%% 

clearvars -except nfish path VSDI VSDmov movie_ref
condi = 404;% ATT: the threshold will be computed respect to the maxim um condition
nsample = 4; %@ SETTTT

plottiles = 1; %also plots early-peak
savetiles = 1;

plot_earlypeak = 0;
save_earlypeak = 0;

% FOT TILES AND EARLY-PEAK FRAMES
fact_thresh =0.4; % @SET : limits parameters
fact_clim= 1.2;


savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/def_figs/tiles/'; % CAHNGEEEEEE

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on= 3;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.forcein = 0; %

%----------------------------------------------------------------
% @SET: for roi-waves plot
%----------------------------------------------------------------

roikind = 'circle'; %
% roikind = 'anat';

%%COMPUTE SETTINGS

%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS
%----------------------------------------------------------------

rejectidx  = compute_rejectidx(VSDI, reject_on, setting);


% GET WAVE AND PLOT 
[sel_trials] = find(VSDI.condition(:,1)==condi);
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
    end
    
    back = VSDI.backgr(:,:,VSDI.nonanidx(1));
    
    %to plot single trial
    movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4);
    movie2plot(:,:,end) = back; %clean non-blured background
   
    
% GET ROI, COORDINATES AND MASK
close all

for i = 1:nsample
newroi_label{i} =num2str(i);
end
% nsample = numel(newroi_label); 
r = VSDI.roi.circle.R;
[coord, roimask] = roicir_draw(VSDI.crop.preview,newroi_label,r); 


% ---------------------------------------------------------------
    % GET WAVES
    %----------------------------------------------------------------
        % -------------------------------------------------------
        % CALCULATE %F WAVE FOR EACH ROI
        % -------------------------------------------------------
        meanF0 = squeeze(mean(VSDmov.F0(:,:,sel_trials),3));
        
            %                 allroi_waves(:,roii,ci) = roi_TSave(movieave,roimask);
            for i = 1:nsample
            waveroi(:,i) = roi_TSave_percF_roiwise(movie2plot,roimask(:,:,i), meanF0);
            end
        % -------------------------------------------------------
        % RANGE
        % -------------------------------------------------------
        
        wave.start_ms = 0; % time in ms for first tile
        wave.end_ms = 1200;
%        wave.start= dsearchn(VSDI.timebase, wave.start_ms);
%         wave.end = dsearchn(VSDI.timebase, wave.end_ms);
        
        
        % PLOT
        figure
        cmap = roicolors();
        cmap = cmap(1:2:end,:); %roicolors map has double values for 2 hemispheres
        
        back = VSDI.backgr(:,:,VSDI.nonanidx(1));
        
        sp1 = subplot(1,2,1);
        switch roikind
            case 'circle'
                roicirc_preview_multiple_cmap(back, coord, VSDI.roi.circle.R, sp1, cmap);
        end
        
        sp1.Visible = 0;
        
        sp2= subplot(1,2,2);
        
        pos1 = get(sp1,'Position');
        pos2 = get(sp2,'Position');
        pos3= [pos2(1) pos2(2) pos1(3) pos1(4)];
        set(sp2, 'Position',pos3)
        
%         nroi = length(selroi);
%         roicolors= roi_colors();
        
        hold on
       
        for i = 1:nsample
            plot(VSDI.timebase, waveroi(:,i), 'linewidth', 1.3, 'Color', cmap(i,:));
        end
        
        % create legend (coord)
        for i = 1:nsample
           L{i} = ['[' num2str(round(coord(i,1))) ',' num2str(round(coord(i,2))) ']'];
        end
        
        legend(L, 'Location', 'northoutside')
sgtitle([num2str(VSDI.ref) '. condi:' num2str(condi)])
