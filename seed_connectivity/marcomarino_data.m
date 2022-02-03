
% ref_movie= '_18filt6';
nfish=11;
ref_movie= '_16diff_f0pre';
% ref_movie = '_18filt6';

VSDI = TORus('load', nfish);
VSDmov = TORus('loadmovie',nfish,ref_movie);


condition = 402;

sel_trials = find(VSDI.condition(:,1)==condition);

% sel_trials = setdiff(sel_trials, rejectidx);

% data.filtered = VSDmov.data(:,:,1:end-1,sel_trials);
% data.backgr =  VSDI.backgr(:,:,VSDI.nonanidx(1));
% data.def = 'n11 - c402 - 16filt 18filt';
data.unfiltered = VSDmov.data(:,:,1:end-1,sel_trials);

movieave = mean(VSDmov.data(:,:,:,sel_trials),4);

act = VSDmov.data(:,:,1:end-1,sel_trials);
maxval = max(act(:));
% tileset.clims = [-maxval*0.6 maxval*0.6];
% tileset.thresh = [-maxval*0.1 maxval*0.1];
tileset.clims = [ ];
tileset.thresh = [0 0];

tileset.custom_map = jet;

tileset.nrowcol = [4 6];
tileset.backgr = VSDI.backgr(:,:,VSDI.nonanidx(1));
% tileset.interp =6;
tileset.start_ms = -20;
tileset.end_ms = 1000;

plot_tilemovie_custom(movieave, VSDI.timebase, tileset, [], jet);

sgtitle('average movie: unfiltered data')


clearvars -except rejectidx data