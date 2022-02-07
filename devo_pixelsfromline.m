% GET PIXELS THAT CROSS A LINE 

user_settings
VSDI = TORus('load', 11) 

im = VSDI.backgr(:,:,VSDI.nonanidx(1));
colormap(bone)
imagesc(im)
axis image 

roi = drawline()

% get coord
xs = [roi.Position(1,1) roi.Position(2,1)];
ys = [roi.Position(1,2) roi.Position(2,2)];

% get pix values 
pixval = improfile(im, xs, ys);

figure
% plot(pixval, 1:length(pixval)); axis tight; title('pixel values')
plot(1:length(pixval),pixval); axis tight; title('pixel values')


% Draw line 

% Get background line F0line

% Get activity dF pixels from 18filt movie Aline

% Convert to %F  Aline./F0line


%% MAKE RASTER OF LINE ALONG TIME

%% MAKE RASTER OF LINE ALONG TRIALS 