% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

clear
load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot.mat')

% ///////////////////////////////////////////////////////////
% SETTINGS

selroinames = {'dm4m_R',  'dm2_R'};%dm3,
roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_17filt5' ;
% ref_movie= '_12filt5' ;

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures' ;%@ SET
ploton = 0;
plotbar_on = 1;
% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1434]; %ms

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 600]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-25 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350];

% feedf.noise.fr_abovenoise = 30;
% feedf.noise.SDfactor = 2;

% Window For average-based analysis

feedf.window_ave = feedf.window;

% feedf.noise_ave = feedf.noise;
% feedf.noise_ave.SDfactor = 4;% SET differences


feedf.method = 'movsum';
% END OF FUNCTION SETTINGS

% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{4,1} = 'wmean (ms)';
params{4,2} = feedf.window.wmean;
% END


% ///////////////////////////////////////////////////////////////

%----------------------------------------------------------------
% @SET: REJECT SETTINGS
%----------------------------------------------------------------

reject_on = 1;  %@ SET
% Subsettings:
setting.manual_reject = 0; %@ SET
setting.GSmethod_reject = 1;  %@ SET
setting.GSabsthres_reject = 1; %@ SET+
setting.force_include = 0; %@ SET


% USER SETTINGS
path.rootpath = '/home/tamara/Documents/MATLAB/VSDI/TORus'; %@ SET for each experiment
tempsep.idcs   = strfind(path.rootpath,'/');
tempsep.newdir = path.rootpath(1:tempsep.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data','BVlists');

addpath(genpath(fullfile(tempsep.newdir, 'VSDI_ourToolbox', 'functions')));
% END USER_SETTINGS

% BUILD MEASURES GROUP -LONGFORMAT- MATRIX
i = 1;
for suji  =  1:size(groupplot,1)

    nfish = groupplot{suji,1};
    cond_lohi = groupplot{suji,3};

    VSDI = TORus('load', nfish);
    VSDmov = TORus('loadmovie',nfish,ref_movie);
    movies = VSDmov.data ;


    %----------------------------------------------------------------
    % CONTROL ROI PICTURE
    %----------------------------------------------------------------
    %     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R);
    %     title([num2str(VSDI.ref) 'roi preview:' selroinames{:}])
    %     saveas(gcf, [num2str(VSDI.ref)'roipreview'], 'jpg')

    %----------------------------------------------------------------
    % SELECT ROI
    %----------------------------------------------------------------

    switch roikind
        case 'circle'
            selroi =name2idx(selroinames, VSDI.roi.labels_circ);
        case 'anat'
            selroi =name2idx(selroinames, VSDI.roi.labels);
    end

    %----------------------------------------------------------------
    % GET INDEXES OF TIMERANGE
    %----------------------------------------------------------------
    idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
    idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values

            %----------------------------------------------------------------
            % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
            %----------------------------------------------------------------

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


            % -------------------------------------------
            % FIRST LOOP THROUGH CONDITIONS AND ROIS TO GET THRESHOLD FOR LATENCY RESPECT TO THE MAX PEAK (of all regions, i.e., a
            % common threshold to all)
                % -------------------------------------------
                % we get and store the value in this first loop to get the max
                % value, and in the second loop we use that max-val as
                % threshold

            thesh = [];
            for condi =   1:3

                condition = cond_lohi(condi);

                % -------------------------------------------
                % SELECT CASES  AND AVERAGE MOVIE
                % -------------------------------------------
                sel_trials = find(VSDI.condition(:,1)==condition);
                if reject_on
                    sel_trials = setdiff(sel_trials, rejectidx);
                end

                movieave = mean(movies(:,:,idxrange,sel_trials),4);
                movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                ...(that normally corresponds to the background)


                for roi_i =  makeRow(selroi)

                    % STORE THE WAVES TO AVOID RECOMPUTE IT IN THE NEXT LOOP

                    switch roikind
                        case 'circle'
                            roimask = VSDI.roi.circle.mask(:,:,roi_i);
                            %             anamask = VSDI.roi.manual_mask(:,:,1);
                        case 'anat'
                            roimask = VSDI.roi.manual_mask(:,:,1);
                    end

                % STORE WAVES. Use timerange set

                    roiwave =  roi_TSave(movieave,roimask);
                    roiwaves(:,roi_i, condi, suji) = roiwave;
%                     testwave = roi_TSave(movieave,roimask);
%                     plot(testwave)

                    % GET THE THRESHOLD to
                    thesh(roi_i, condi)= max(movmean(roiwave,10)); %20% of the maximum value


                end % second loop of roi_i

               end % condi


                % -------------------------------------------------------------------
                % -------------------------------------------------------------------

                % --------------------------------------------------------------------
                % SECOND LOOP THROUGH CONDI-ROIS TO GET THE MEASURE OF ABSOLUTE THRESHOLD AND STORE THEM IN LONG FORMAT
                % --------------------------------------------------------------------

            for condi = 1:3

                condition = cond_lohi(condi);

                % -------------------------------------------
                % SELECT CASES  AND AVERAGE MOVIE
                % -------------------------------------------
                sel_trials = find(VSDI.condition(:,1)==condition);
                if reject_on
                    sel_trials = setdiff(sel_trials, rejectidx);
                end

                movieave = mean(movies(:,:,idxrange,sel_trials),4);
                movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                ...(that normally corresponds to the background)

                for roi_i = makeRow(selroi)
                    % -------------------------------------------
                    % GET ROI WAVES and MASKS AGAIN
                    % -------------------------------------------

                    switch roikind
                        case 'circle'
                            roimask = VSDI.roi.circle.mask(:,:,roi_i);
                        case 'anat'
                            roimask = VSDI.roi.manual_mask(:,:,1);
                    end

                    roiwave =roiwaves(:,roi_i, condi, suji); % GET BACK FROM PREVIOUS LOOP

                    % -------------------------------------------
                    % %dF/F0
                    % -------------------------------------------

                    % plot(roiwave(1:601))

                    % -------------------------------------------
                    % CALCULATE MEASURES
                    % -------------------------------------------

                    % peakminubaseline
                    outputP = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window,[], feedf.method, 0, 0);
%                     plot(VSDI.timebase(idxrange),roiwave)
                    % initialmean

                    % latency respect to the 30% of the own peak

                    % -------------------------------------------
                    % LATENCY 30% ONSET
                    % -------------------------------------------
                    % find the threshold

                    output_lat = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window,[], feedf.method, 0, 0);


                    % -------------------------------------------
                    % LATENCY RESPECT  THE 20% OF ALL PEAKS
                    % -------------------------------------------
                    % to find the risingtime over the absolute threshold
                    maxvalue = max(thesh(:));
                    feedf.window.risingthresh = maxvalue*0.2;
                    % specific inputs to use an absolute threshold
                    output_abs = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window,[], feedf.method, 0, 0);


                    % ibi0

                    % -------------------------------------------
                    % STORE IN LONG FORMAT
                    % -------------------------------------------

                    % FACTORS
                    longF{i,1} = suji ; %subject id
                    longF{i,2} = roi_i; % roi
                    longF{i,3} =VSDI.roi.labels;
                    longF{i,4} = condi; % condition (1 = low; 2 = hi)
                    % MEASURES
                    longF{i,5} = round(outputP.peakminusbasel,2); %
                    longF{i,6} = outputP.peaklatency; %
                    longF{i,7} = output_lat.onset30_latency_ms; %
                    longF{i,8} = output_abs.risingthreshlat; %
                    longF{i,9} = round(outputP.wmean,2);%
                    
                    % -------------------------------------------
                    % ... AND NORMALIZED
                    % -------------------------------------------
                    
                    longFnorm(i,1) = suji ; %subject id
                    longFnorm(i,2) = roi_i; % roi
                    longFnorm(i,3) = condi; % condition (1 = low; 2 = hi)
                    % MEASURES
                    longFnorm(i,4) = round(outputP.peakminusbasel,2); %
                    longFnorm(i,5) = outputP.peaklatency; %
                    longFnorm(i,6) = output_lat.onset30_latency_ms; %
                    longFnorm(i,7) = output_abs.risingthreshlat; %
                    longFnorm(i,8) = round(outputP.wmean,2);%

                    
                    
                    
                    i = i+1;

                    % -------------------------------------------
                    % STORE IN WIDE FORMAT (for grouped barplots)
                    % -------------------------------------------
                    wideF.peakminusbasel(roi_i, condi) = round(outputP.peakminusbasel,2);
                    wideF.peaklatency(roi_i, condi)  =  outputP.peaklatency; %
                    wideF.onset30_latency_ms(roi_i, condi) = output_lat.onset30_latency_ms; %
                    wideF.risingthreshlat(roi_i, condi) = output_abs.risingthreshlat; %
                    wideF.wmean(roi_i, condi) = round(outputP.wmean,2);%

                    clear outputP output_lat output_abs

                end % roi_i
                clear  roiwave thresh movieave  sel_trials


            end % condi

            clear reject_idx

                % --------------------------------------------------------
                % PLOTS : HARDWIRED FOR 2-ROIS AND 2-CONDITIONS !!! plot both rois from both conditions (in each
                % subject)
                % --------------------------------------------------------
                if ploton
                figure
                ax(1) = subplot(2,3,[1 2 4 5]);
                plot(VSDI.timebase(idxrange), roiwaves(:,selroi(1), 1, suji), '--b', 'linewidth', 1, 'DisplayName',strcat(selroinames{1}, '-low'))
                hold on
                               plot(VSDI.timebase(idxrange), roiwaves(:,selroi(1), 2, suji), '-b', 'linewidth', 1, 'DisplayName',strcat(selroinames{1}, '-high'))
                               plot(VSDI.timebase(idxrange), roiwaves(:,selroi(1), 3, suji), ':b', 'linewidth', 1, 'DisplayName',strcat(selroinames{1}, '-control'))

                                plot(VSDI.timebase(idxrange), roiwaves(:,selroi(2), 1, suji), '--r',  'linewidth', 1, 'DisplayName',strcat(selroinames{2}, '-low'))
                               plot(VSDI.timebase(idxrange), roiwaves(:,selroi(2), 2, suji), '-r',  'linewidth', 1, 'DisplayName',strcat(selroinames{2}, '-high'))
                               plot(VSDI.timebase(idxrange), roiwaves(:,selroi(2), 3, suji), ':r', 'linewidth', 1, 'DisplayName',strcat(selroinames{2}, '-control'))


                xlim([VSDI.timebase(idxrange(1)), VSDI.timebase(idxrange(end))])
                legend()
                tname = [num2str(VSDI.ref) '(' ref_movie(2:end) '-' roikind ')' ];

                ax(2) = subplot(2,3,3);

                    switch roikind
                        case 'circle'
                roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R, ax(2));
                        case 'anat'
                           roi_preview_multiple(VSDI.crop.preview,  VSDI.roi.manual_poly(selroi,:), ax(2))
                    end
                    ax(2).Visible = 'off';
                %     title([num2str(VSDI.ref) 'roi preview:' selroinames{:}])


                sgtitle(tname)
                saveas(gcf, fullfile(savein, ['plot' tname  '_reject' num2str(reject_on) '.jpg']), 'jpg')
                close
%                 --------------------------------------------------------
                end %ploton

                % --------------------------------------------------------
                % BARPLOTS : HARDWIRED FOR 2-ROIS AND 2-CONDITIONS !!! plot both rois from both conditions (in each
                % subject)
                % --------------------------------------------------------
                if plotbar_on

                   sgtitle(['fish-' num2str(VSDI.ref)])
                   subplot(2,3,1)
                   bar(wideF.peakminusbasel(selroi,:))
                   title('peak')

                   subplot(2,3,2)
                   bar(wideF.peaklatency(selroi, :)); %
                   title('peak latency')

                   subplot(2,3,3)
                   bar( wideF.onset30_latency_ms(selroi, :)); %
                   title('onset latency 30%')

                   subplot(2,3,4)
                   bar(wideF.risingthreshlat(selroi, :)); %
                   title('onset latency 20%(abs)')

                   subplot(2,3,5)
                   bar(wideF.wmean(selroi, :));%
                   title('wmean (onset lat)')

                   saveas(gcf, fullfile(savein, ['barplots' num2str(VSDI.ref) 'reject' num2str(reject_on) '.jpg']), 'jpg')
                   close
                end

            clear roiwaves

    end % groupplot

    %% -------------------------------------------
    % EXPORT FOR R
    % -------------------------------------------
    excelname = fullfile(savein, ['long_format_forR_group.xls']);
    labels = {'id' 'roi' 'cond' 'peak' 'peaklat' 'lat30own_ms' 'lat20abs_ms', 'wmean'};

    % write output (new sheet for each fish
    writematrix (longF, excelname, 'sheet',  [roikind ref_movie 'rej' num2str(reject_on)])
    writecell (labels, excelname, 'sheet', 'labels')
    writecell (groupplot, excelname, 'sheet','group')
    writecell (params, excelname, 'sheet','param')

    
    
% %% BARPLOT
% labelcondi = {'low', 'high'};
% for suji =3% 1:length(groupplot)
%     for condi= 1:2
%     for roi_i = 1:7
%         datapeak(condi,roi_i) =  bar.peak(:,roi_i, condi, suji) ;
%     end
%     ax(condi) = subplot(1,2,condi);
%     bar(ax(condi), datapeak);
%     title(labelcondi(condi));
%     end
% end

blob()
