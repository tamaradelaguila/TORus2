clear
user_settings

%----------------------------------------------------------------
% SETTINGS
%----------------------------------------------------------------
nfish = 11;
conditions = [400:404];
ref_movie = '_16diff_f0pre';

exclude_basel = 0;
% ncomp = 15;

pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/ICA/';
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
% TRANFORM TO NIFTI
%----------------------------------------------------------------
filename = fullfile(pathsave, 'NIFTItest');
niftiwrite(supermovie,filename)

% V = niftiread(filename);
info = niftiinfo(filename);

V = supermovie(:,:,1:100); 
niftiwrite(V,filename)
