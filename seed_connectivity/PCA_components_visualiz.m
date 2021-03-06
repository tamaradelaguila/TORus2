clear
user_settings

%----------------------------------------------------------------
% SETTINGS
%----------------------------------------------------------------
nfish = 11;
conditions = [1000 1002 1003];
ref_movie = '_16diff_f0pre';

exclude_basel = 0;
ncomp = 15;

pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/seed_connectivity/PCA_removal';
%----------------------------------------------------------------
% LOAD DATA
%----------------------------------------------------------------
VSDI = TORus('load', nfish);
VSDmov= TORus('loadmovie', nfish, ref_movie);

idxcondi = [];
for ci = 1:numel(conditions)
    temp = choose_condidx(VSDI.condition(:,1),conditions(ci));
    idxcondi = [idxcondi temp];
end
idxcondi = makeCol(sort(idxcondi));

%----------------------------------------------------------------
% %  CONCATENATE
%----------------------------------------------------------------

if exclude_basel
    t0 = dsearchn(VSDI.timebase, 0);
    moviesin= VSDmov.data(:,:,t0:end-1,idxcondi); %from t0 to the end (excluding background)
else
    moviesin= VSDmov.data(:,:,1:end-1,idxcondi); %from t0 to the end (excluding background)
end

dim = size(moviesin);
newdim = dim(3)*dim(4);

supermovie = reshape(moviesin, [dim(1) dim(2) newdim]); %concatenated

%----------------------------------------------------------------
% PCA
%----------------------------------------------------------------

dim2 = size(supermovie);
X = reshape(supermovie, [dim2(1)*dim2(2) dim2(3)]); %pixels in a row
X = X'; % observations in rows

%De-mean:
mu = mean(X);
X = bsxfun(@minus,X,mu);

% PCA:
[eigenvectors, scores] = pca(X);

%----------------------------------------------------------------
% ICA
%----------------------------------------------------------------

dim2 = size(supermovie);
X = reshape(supermovie, [dim2(1)*dim2(2) dim2(3)]); %pixels in a row
X = X'; % observations in rows
% 
% %De-mean:
% mu = mean(X);
% X = bsxfun(@minus,X,mu);

% PCA:
q = 15; % nº of feats to extract
Mdl = rica(X,q); 

z = transform(Mdl,X); %transforms the data x into the features z via the model Mdl.

    Rmovies =reshape(z, dim2);
    Rmovies = reshape(Rmovies, dim);

%% ----------------------------------------------------------------
% COMPONENTS RECONSTRUCTION and PLOT TILES
%----------------------------------------------------------------

for n = 1:ncomp
    %  RECONSTRUCT
    % Data reconstructed from certain components
    Xhat = scores(:,n) * eigenvectors(:,n)';
    Xhat = bsxfun(@plus, Xhat, mu);
    Rmovies =reshape(Xhat', dim2);
    Rmovies = reshape(Rmovies, dim);
    
    % TILE OF COMPONENT-RECONSTRUCTED MOVIE (condition-specific)
    for condi = makeRow(conditions)
        idx = choose_condidx(VSDI.condition(:,1),condi); 
        
        % get idx relative to those used to compute the PCA
        idxrel = dsearchn(idxcondi,makeCol(idx));%it has to be relative to the computed trials
        
        % get movies from the condition and average
        movie = mean(Rmovies(:,:,:,idxrel),4);
        
        tileset.start_ms = -10;
        tileset.end_ms = 800;
        
        maxval = max(movie(:));
        tileset.thresh = [0 0];
        tileset.clims = [];
        tileset.nrowcol = [1 6];
        tiles.backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));
        
        plot_tilemovie_custom(movie, VSDI.timebase, tileset);
        tit= [num2str(VSDI.ref) ':PC' num2str(n) ' - cond' num2str(condi) ];
        sgtitle(tit)
        
        saveas(gcf, fullfile(pathsave, [tit '.jpg']), 'jpg') % prints it as you see them
        close
        idx = []; idxrel = [];
    end % for condi
    
end % for n