%% SEED-BASED connectivity

% The PCA will be performed on all EI + blank trials (leaving out the tone
% conditions) and the trials from each condition will be reconstructed from
% all-but-the-first component

cleardata.filtered
user_settings
%----------------------------------------------------------------
% @SET: BASIC PARAMETERS
%----------------------------------------------------------------

savein = '/home/tamara/Documents/MATLAB/VSDI/TORus/seed_connectivity/seedCONN_1PCremoval' ;%@ SET
% load('/home/tamara/Documents/MATLAB/VSDI/STIM/PCAgroupinfo.mat') %stores group parameters
% load('/home/tamara/Documents/MATLAB/VSDI/TORus/plot/informes/03_figure_sketch/groupplot_measures/groupplot_new.mat') % ATTENTION: select cases manually (there are repetitions)

seednames = {'dm4m_R2', 'dm4_R', 'dm2_R2', 'dm1_R', 'dldm_R2', 'dm3_R2', 'dldr_R',...
            'dm4m_L2', 'dm4_L', 'dm2_L2', 'dm1_L', 'dldm_L2', 'dm3_L2', 'dldr_L',};

roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_16diff_f0pre';% '_17filt5' ; '_18filt6'; '_19filt7'; '_16diff_f0pre'

activ_unit = 'diffF'; % @ MANUALLY SET (just for display purposes)
analysisref = 'new_group5_'; % MANUALLY SET!!! extra info for the name. Set group according to the rows of 'groupplot' selected in;   for suji  =  [1 3 4 9]


%----------------------------------------------------------------
%  @SET: SUBJECT-SPECIFIC PARAMETERS
%----------------------------------------------------------------

nfish = 12;
conditions = 400:404;
tiles_on = 1;

%----------------------------------------------------------------
% @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
%----------------------------------------------------------------
reject_on= 3;

setting.manual_reject = 1; %
setting.GSmethod_reject = 1;  %
setting.GSabsthres_reject = 1; %
setting.force_include = 0; %



%% PERFORM PCA ON ALL TRIALS FROM ALL CONDITIONS

% In each loop, the parameters will be extracted from the structure.
% In subsequents loops, the seed-correlations, PCA and plots will be
% performed

%----------------------------------------------------------------
% GET PARAMETERS FROM STRUCTURE
%----------------------------------------------------------------

exclude_basel = 0;

clean_PCA = 0;
clean_seed = 1; % can be zero only with: clean_PCA = 0;

%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus('load', nfish);
VSDmov = TORus('loadmovie', nfish, ref_movie);


%----------------------------------------------------------------
% COMPUTE REJECTION IDX FROM REJECT-OPTIONS 
%----------------------------------------------------------------
rej = 'reject' ;
if reject_on > 1
rej = [rej num2str(reject_on)];
end

rejectidx = [];


if setting.manual_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
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

%----------------------------------------------------------------
% IDX SELECTION 
%----------------------------------------------------------------

idxsel = [];
for ci = 1:numel(conditions)
    temp = choose_condidx(VSDI.condition(:,1),conditions(ci));
    idxsel = [idxsel temp];
end


if clean_PCA
  idxsel = setdiff(idxsel, rejectidx);   
end

idxsel = sort(idxsel');

%----------------------------------------------------------------
% GET MOVIE DATA - now relative to the indexes 'idxsel'
%----------------------------------------------------------------

if exclude_basel
    t0 = dsearchn(VSDI.timebase, 0);
    movies_in= VSDmov.data(:,:,t0:end-1,idxsel); %from t0 to the end (excluding background)
else
    movies_in= VSDmov.data(:,:,1:end-1,idxsel); %from t0 to the end (excluding background)
end


dim = size(movies_in);
newdim = dim(3)*dim(4);

supermovie = reshape(movies_in, [dim(1) dim(2) newdim]); %concatenated


%----------------------------------------------------------------
% PCA: 1PC SUBSTRACTION
%----------------------------------------------------------------

        dim2 = size(supermovie);
        X = reshape(supermovie, [dim2(1)*dim2(2) dim2(3)]); %pixels in a row
        X = X'; % observations in rows
        
        %De-mean:
        mu = mean(X);
        X = bsxfun(@minus,X,mu);
        
        % PCA:
        [eigenvectors, scores] = pca(X);
        
        % Data reconstructed from certain components
        Xhat_butpc1 = scores(:,2:15) * eigenvectors(:,2:15)';
        Xhat_butpc1 = bsxfun(@plus, Xhat_butpc1, mu);
        new_supermovie =reshape(Xhat_butpc1', dim2);

        movies_new = reshape(new_supermovie, dim);
        

%% FROM EACH CONDITION: SUBSTRACT AND COMPUTE SEED-BASED CONNECTIVITY
j = 0;
for condi = conditions
j = j+1;

% Get index relative to the already selected movies 
idxcond = choose_condidx(VSDI.condition(:,1),condi);
if clean_seed
  idxcond = setdiff(idxcond, rejectidx);   
end

idxnew = dsearchn(idxsel,idxcond'); %relative to the reconstructed movies (movies_new) but also to movies_in

% Get movies from the selected condition
movie = movies_in(:,:,:,idxnew); 
dim3 = size(movie); 
movieconc = reshape(movie, [dim3(1) dim3(2) dim3(3)*dim3(4)]) ;

movie_new = movies_new(:,:,:,idxnew); 
dim3 = size(movie_new); 
movieconc_new = reshape(movie_new, [dim3(1) dim3(2) dim3(3)*dim3(4)]) ;

cond_labels{j,1} = condi;
i = find(VSDI.condition(:,1) == condi, 1, 'first');
cond_labels{j,2} = [num2str(VSDI.condition(i,4)) 'mA'];
clear i

%----------------------------------------------------------------
% SEED-BASED CORRELATION OF ORIGINAL MOVIE
%----------------------------------------------------------------
seed_allidx = name2idx(makeCol(seednames), VSDI.roi.labels_circ);

for seedi = makeRow(seed_allidx)
    seedname = VSDI.roi.labels_circ{seedi};
    mask =  VSDI.roi.circle.mask(:,:,seedi);
    %         centre = seeds.centre(seedi,:);
    %         r = seeds.r;
    
    seedwave = roi_TSave(movieconc,mask);
    seedwave = seedwave';
    
    for rowi =1:dim(2)
        for coli = 1:dim(1)
            pixelwave = squeeze(movieconc(coli, rowi, 1:end-1)); %roi_TSave substracts the last frame so we have to match length
            
            if sum(pixelwave) ==0
                rhomap(coli, rowi, j, seedi) = 0;
            else
                rhomap(coli, rowi, j, seedi) = corr(seedwave, pixelwave);
            end
            
            clear pixelwave
        end
    end % rowi
clear seedwave
end % for seedi

% 
% %----------------------------------------------------------------
% % % PLOT TILES: ORIGINAL DATA
% %----------------------------------------------------------------
% tileset.start_ms = 0; % time in ms for first tile
% tileset.end_ms = 800;
% tileset.time2plot = 0; %select time (ms)
% tileset.x = 35;
% tileset.y = 35;
% tileset.nrowcol = [6 4];
% 
% if tiles_on
%     
%     ti = ['(original) cond:' num2str(cond(condrow))];
% end
% 
%     % -------------------------------------------
%     % get inputs to the function
%     movieplot = mean(movies_in,4);
%     savetiles{j} = movieplot; %to further inspect
%     
%     maxval = max(movieplot(:));
%     maxval = maxval;
%     tileset.clims = [-maxval maxval] ;
%     tileset.thresh = [-maxval/6 maxval/6];
%     tileset.backgr = VSDI.backgr(:,:,idxcondi(1));
%     
%     % -------------------------------------------
%     plot_tilemovie_custom(movieplot, VSDI.timebase, tileset,[])
%     
%     sgtitle(ti)
%     name = [ 'n' num2str(groupRow) '-' experiment '-' num2str(VSDI.ref) ref_movie 'TILES ' num2str(j) '-A-original.jpg'];
%     name2save = fullfile(savein, name);
%     saveas(gcf, name2save, 'jpg')
%     close
%     

%----------------------------------------------------------------
% SEED-BASED CORRELATION OF RECONSTRUCTED DATA: PC-1 out
%----------------------------------------------------------------
for seedi = makeRow(seed_allidx)
    seedname = VSDI.roi.labels_circ{seedi};
    mask =  VSDI.roi.circle.mask(:,:,seedi);
    %         centre = seeds.centre(seedi,:);
    %         r = seeds.r;
    
    seedwave = roi_TSave(movieconc_new,mask);
    seedwave = seedwave';
    
    for rowi =1:dim(2)
        for coli = 1:dim(1)
            pixelwave = squeeze(movieconc_new(coli, rowi, 1:end-1)); %roi_TSave substracts the last frame so we have to match length
            
            if sum(pixelwave) ==0
                rhomap_new(coli, rowi, j, seedi) = 0;
            else
                rhomap_new(coli, rowi, j, seedi) = corr(seedwave, pixelwave);
            end
            
            clear pixelwave
        end
    end % rowi

clear seedwave
end % for seedi

% %----------------------------------------------------------------
% % SPATIAL SMOOTH OF SEED-BASED CORRELATION MAPS
% %----------------------------------------------------------------
% 
% pix = 3; %smoothing kernel
% 
% for seedi = 1:size(seeds.centre,1)
%     maps = squeeze(rhomap_butpc1(:,:,:,seedi));
%     
%     kernel = ones(pix);
%     
%     for fr = 1:4
%         frame = squeeze(maps(:,:,fr));
%         kernel = kernel / sum(kernel(:));
%         localAverageImage = conv2(double(frame), kernel, 'same');
%         
%         rhomap_butpc1_smoothed(:,:,fr,seedi) = localAverageImage;
%     end
%     clear maps
% end


% 
% % % PLOT TILES: OF RECONSTRUCTED DATA: EXCEPT PC1
% %----------------------------------------------------------------
% if tiles_on
%     
%     % get condition info
%     switch experiment
%         case 'TORus'
%             ti = ['(noPC1) cond: ' num2str(cond(condrow))];
%         case 'STIM'
%             ti = ['(noPC1) p' num2str(VSDI.condition(idxcondi(1),1)) '-' num2str(VSDI.condition(idxcondi(1),2)) '\mum -' num2str(VSDI.condition(idxcondi(1),3)) '\muA'];
%         case 'TORtone'
%             ti = ['(noPC1) hz: ' num2str(cond(condrow))];
%             
%     end
%     % -------------------------------------------
%     % get inputs to the function
%     movieplot = mean(movie_butpc1,4);
%     savetiles_butpc1{j} = movieplot; %to further inspect
%     maxval = max(movieplot(:));
%     maxval = maxval;
%     tileset.clims = [-maxval maxval] ;
%     tileset.thresh = [-maxval/6 maxval/6];
%     tileset.backgr = VSDI.backgr(:,:,idxcondi(1));
%     
%     % -------------------------------------------
%     plot_tilemovie_custom(movieplot, VSDI.timebase, tileset,[])
%     
%     sgtitle(ti)
%     name = [ 'n' num2str(groupRow) '-' experiment '-' num2str(VSDI.ref) ref_movie 'TILES ' num2str(j) '-C-noPC1 reconstr.jpg'];
%     name2save = fullfile(savein, name);
%     saveas(gcf, name2save, 'jpg')
%     close
%     
% end %if tiles_on


clear movie movie_new movieconc moviecond_new idxcond idxnew
end
clear j


labels{1} = 'pix';
labels{2} = 'pix';
labels{3} = 'cond';
labels{4} = 'seeds';

seedCONN.ref = VSDI.ref;
seedCONN.rho = rhomap;
seedCONN.rho_no1pc = rhomap_new ; 
seedCONN.labels = labels;
seedCONN.cond = cond_labels ;
seedCONN.seednames = seednames;

matname = ['RHOmaps_PCA' num2str(clean_PCA),'_seed' num2str(clean_seed), '_' num2str(VSDI.ref) ref_movie];
save(fullfile(savein, matname), 'seedCONN')
blob()

%% PLOT ALREADY COMPUTED SEED-BASED CONNECTIVITY RHOMAPS
% For each seed, plot all conditions and

% EXTRACT DATA AND INFO FROM SAVED RHO-CONNECTIVITY MAPS
% matname = ['RHOmaps' num2str(VSDI.ref) ref_movie]; %simple name (old)
matname = ['RHOmaps_PCA' num2str(clean_PCA),'_seed' num2str(clean_seed), '_' num2str(VSDI.ref) ref_movie];
load(fullfile(savein, matname))
rhomap = seedCONN.rho;
rhomap_new = seedCONN.rho_no1pc ; 
labels = seedCONN.labels ;
cond_labels =seedCONN.cond ;
seednames = seedCONN.seednames;


seed_allidx = name2idx(makeCol(seednames), VSDI.roi.labels_circ);

for seedi = seed_allidx
    name = VSDI.roi.labels_circ{seedi};
    centre = VSDI.roi.circle.center(seedi,:);
    r =  VSDI.roi.circle.R;
    
    %----------------------------------------------------------------
    % % % #2a: PLOT RHO-MAP: ORIGINAL DATA
    %----------------------------------------------------------------
    
    figure
    for ploti = 1:length(cond_labels)

        % ------------------------------
        subplot(1,5,ploti);
        imagesc(rhomap(:,:, ploti, seedi))
        colormap(flipud(hot))
        axis image
        set(gca, 'clim', [0 1])
        title(cond_labels{ploti,2})
        colorbar
    end % ploti
    
    % PLOT SEED ROI IN THE LAST PLOT
    
    ax1= subplot(1,5,5);
    backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));
    imagesc(backgr); colorbar; hold on
    colormap(ax1, 'bone')
    axis image
    title('seed:')
    viscircles(centre, r, 'color','r')
    axis image
    set(gcf, 'Position', get(0, 'Screensize'));
    
    
    supertit = [num2str(VSDI.ref) ref_movie 'Seed-connectivity:' num2str(seedi) '(' name ')'  '-A-original.jpg'];
    sgtitle(supertit)
%     name2save = fullfile(savein, supertit);
    name2save = fullfile(savein, ['PCA' num2str(clean_PCA) '_seed' num2str(clean_seed) '_' supertit]);
    saveas(gcf, name2save, 'jpg')
    close
    
    %----------------------------------------------------------------
    % % % #2c: PLOT RHO-MAP: PC1 OUT
    %----------------------------------------------------------------
    
  figure
    for ploti = 1:length(cond_labels)

        % ------------------------------
        subplot(1,5,ploti);
        imagesc(rhomap_new(:,:, ploti, seedi))
        colormap(flipud(hot))
        axis image
        set(gca, 'clim', [0 1])
        title(cond_labels{ploti,2})
        colorbar
    end % ploti
    
    % PLOT SEED ROI IN THE LAST PLOT
    
    ax1= subplot(1,5,5);
    backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));
    imagesc(backgr); colorbar; hold on
    colormap(ax1, 'bone')
    axis image
    title('seed:')
    viscircles(centre, r, 'color','r')
    axis image
    set(gcf, 'Position', get(0, 'Screensize'));
    
    supertit = [num2str(VSDI.ref) ref_movie 'Seed-connectivity:' num2str(seedi) '(' name ')'  '-B-pc1 substracted.jpg'];
    sgtitle(supertit)
%     name2save = fullfile(savein, supertit);
    name2save = fullfile(savein, ['PCA' num2str(clean_PCA) '_seed' num2str(clean_seed) '_' supertit]);
    saveas(gcf, name2save, 'jpg')
    close
    
end %seedi


% %% PLAY AROUND WITH TILES FROM ORIGINAL AND RECONSTRUCTED MOVIES
%             movieplot = savetiles_butpc1{2}  ; %to further inspect
%             maxval = max(movieplot(:));
%             maxval = 20;
%             tileset.clims = [-maxval maxval] ;
%             tileset.thresh = [-maxval/6 maxval/6];
%             tileset.backgr = VSDI.backgr(:,:,idxcondi(1));
%
%             % -------------------------------------------
%             plot_tilemovie_custom(movieplot, VSDI.timebase, tileset,[])
%             sgtitle('pc1-substract reconstructed movie')
%
%             movieplot = savetiles{2}  ; %to further inspect
%             maxval = max(movieplot(:));
%             maxval = 10;
%             tileset.clims = [-maxval maxval] ;
%             tileset.thresh = [-maxval/6 maxval/6];
%             tileset.backgr = VSDI.backgr(:,:,idxcondi(1));
%
%             % -------------------------------------------
%             plot_tilemovie_custom(movieplot, VSDI.timebase, tileset,[])
%             sgtitle('original movie')
%% Created: 25/01/22
% last updated: