function [] = devo_plot_allframes(moviedata,tileset)


% [] = plot_tilemovie(moviedata, times, tileset)
% EXTRACT 12 TILES FROM INPUT MOVIE 3D and STRUCTURE THAT CONTAINS ALL THE SETTINGS

% INPUTS:
%       moviedata: 3D already-inverted data
%       times: timebase corresponding to the different frames (in ms)
%       tileset: structure with settings for tiling in fields: 
%           .start_ms = -300; % time in ms for first tile       
%           .end_ms = 2520;
%           .time2plot = 0; %select time (ms)
%           .x = 35; 
%           .y = 35; 
%           .clims = [] --- values above or below will be saturated
%           .thresh = [negthresh posthres]; --- values inside this range
%         .nrowcol = [nrows ncols] number of rows and columns
%         .step = 10 %frame step in tiles

%           will be zero-out
%           .backgr
%       custom_map 
%       roipoly: roi polygons if we want to overdraw them in every frame (e.g.:VSDI.roi.manual_poly)

% Input control
if ~exist('roipoly', 'var')
    plot_rois = 0;
else
    plot_rois = 1;
end

if nargin < 5
    custom_map = colormap_loadBV();
end

if isfield(tileset, 'backgr')
    background = tileset.backgr;
else
    background = moviedata(:,:,end);
end

if isfield(tileset, 'nrowcol')
nrows = tileset.nrowcol(1) ;
ncols = tileset.nrowcol(2);

else
dipl('the field nrowcol needs to be specified')
end


if isfield(tileset, 'step')
step =tileset.step;
else
step = NaN;
end

% End of input control

% datatime = times; 
% timeindx = 91:340;
imtiles = moviedata(:,:, timeindx); % movie of the selected frames

imback = repmat(3,background);

%% Get alphachannel to plot overlaid

imalpha = imtiles> tileset.thresh(1) & imtiles< tileset.thresh(2);
imalpha = ~imalpha;

%% RESHAPE
 
dim = size(imtiles);
Ncols = dim(3)/nrows ; %number of columns in the final plot is equal to the timeidx divided by the 

flat = [];
flatAlfa=[];
flatB = [];

idx = 1; 
for ii = 1:nrows
    range=idx:idx+Ncols-1;
    
    % ACTIVITY MATRIX
    temp= imtiles(:,:,range);
    s = size(temp);
    temp = reshape(temp, [s(1),s(2)*s(3)]);     % create strip of tiles,
    flat = cat(1, flat,temp); % add to preexisting stripe
    
    
        % alpha
    tempAlfa= imtiles(:,:,range);
    s = size(tempAlfa);
    tempAlfa = reshape(tempAlfa, [s(1),s(2)*s(3)]);     % create strip of tiles,
    flatAlfa = cat(1, flatAlfa,tempAlfa); % add to preexisting stripe

    % BACKGROUND
    tempB = repmat(background, 1, Ncols);
    flatB = cat(1, flatB,tempB); % add to preexisting stripe
    % COUNTER
    idx = Ncols+1;
    
    
    
end

%% PLOT 
figure

   ax1 = subplot(1, 1, 1); 
%    imagesc(imtiles(:,:,ploti))
   plot_framesoverlaid(flat,flatB, flatAlfa ,0, ax1, tileset.clims, 0, tileset.thresh, custom_map); 
   
  
end


%% SUPPORT FUNCTION

function  plot_framesoverlaid(imAct, imBack, logicalpha, plotnow, axH, act_clim, plot_cbar, thresh, custom_map)
% INPUT 
% 'imAct' - image to display in colors
% 'imBack' - background image
% 'logicalpha': logic matrix with transparency information (1 = shown pixel)
% 'plotnow': whether we want to plot the figure (set to 0 if we are providing a axisHandle)
% 'axH': axes handle in which to plot it, when it's going to be included into a subplot
% 'act_clim': color limits for 'imAct'. If no color limits is provided for the activity (color) map, use the max
% and min from the input movie
% 'plot_cbar' : whether we want to plot the colorbar (useful in subplots)


if ~exist('act_clim')
    act_clim = [min(nonzeros(imAct(:))) max(nonzeros(imAct(:)))];
elseif isempty(act_clim)
    act_clim = [min(nonzeros(imAct(:))) max(nonzeros(imAct(:)))];
end

% if ~exist('thresh')
%     thresh = 0;
% elseif isempty(thresh)
%         thresh = 0;
% end

if ~exist('plotnow')
        plotnow= 1;
elseif isempty(plotnow)
        plotnow= 1;
end 

if ~exist('plot_cbar')
        plot_cbar= 1;
elseif isempty(plot_cbar)
        plot_cbar= 1;
end 


if ~exist('custom_map')
        colormode= 0;
elseif isempty(custom_map)
        colormode= 0;
else 
    colormode = 1;
end

% end of input control ------------------------
if plotnow
    fig1 = figure; 
end

thresh_polarmap = min(abs(thresh));

ax1 = axes;
imagesc(imBack);
colormap(ax1,'bone');
ax1.Visible = 'off';
axis image

ax2 = axes;
% imagesc(ax2,imAct,'alphadata',imAct>thresh);
imagesc(ax2,imAct,'alphadata',logicalpha);

    if colormode == 1
    colormap(ax2, custom_map);
    else
    colormap(ax2,polarmap());
    end
    
caxis(ax2, act_clim);
ax2.Visible = 'off';
if plot_cbar
colorbar;
end
axis image

linkprop([axH ax1 ax2],'Position');

% if plotnow
%    close (fig1)
% end

end

%% Created  26/10/21
