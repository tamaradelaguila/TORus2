%% s03 DEFINE ROIS
% Draw ROIs from brain areas.
clear
user_settings
nfish = 8;
VSDI = TORus('load',nfish);

% Draw ROIs and store in structure
VSDI.roi.labels = {'dm4_R', 'dm4_L',...
   'dm4m_R', 'dm4m_L',...
   'dm3_R', 'dm3_L',...
   'dm2_R', 'dm2_L',...
   'dldc_R', 'dldc_L',...
   'dldm_R', 'dldm_L',...
   'dm1_R', 'dmd1_L'}';

VSDI.roi.labels = {'dm4_R',... FOR MAPS EXPERIMENT
   'dm4m_R'...
   'dm2_R'}';

% DEFINE with the help of the function:
[manual_poly, manual_mask] = roi_draw(VSDI.crop.preview,VSDI.roi.labels); %in each column, one set of rois
[manual_poly, manual_mask] = roi_draw(VSDI.backgr(:,:,VSDI.nonanidx(1)),VSDI.roi.labels); %in each column, one set of rois

% View the result
 roi_preview_multiple(VSDI.crop.preview, manual_poly); 
 roi_preview_multiple(VSDI.backgr(:,:,end), manual_poly); %see in a selected frame

% TO REDO A SINGLE ROI, DRAW THEM ALL AND COPY MANUALLY THE COORDINATES
... AND MASK
        [redo_poly, redo_mask] = roi_draw(VSDI.crop.preview,VSDI.roi.labels); 
         roi_preview_multiple(VSDI.backgr(:,:,140), redo_poly); %see in a selected frame
        manual_poly{15,1} = redo_poly{1,1}
        manual_mask(:,:,15:16) = redo_mask(:,:,1:2);
        roi_preview_multiple(VSDI.backgr(:,:,144), manual_poly); %see in a selected frame

VSDI.roi.manual_poly  = manual_poly;
VSDI.roi.manual_mask = manual_mask;

for ii = 1:size(VSDI.roi.labels,1)
    roiname  =VSDI.roi.labels{1,ii};
%     VSDI.roi.labels{2,ii} = roiname_ipsicontra(roiname, VSDI.info.Sside); 
end

TORus('save',VSDI);

roi_preview_multiple(VSDI.crop.preview, VSDI.roi.manual_poly); 

%% STORE IN WAVES-STRUCTURE that we create here (for extraction in '_extract_ROItimeseries):

VSDroiTS.ref = VSDI.ref; 
VSDroiTS.roi = VSDI.roi; 
TORus('savewave', VSDroiTS);

%% CIRCULAR ROIS
% R = 10;
% [manual_poly, manual_mask] = roicir_draw(VSDI.crop.preview,VSDI.roi.labels, R); %in each column, one set of rois
% 
% figure
% imagesc(VSDI.crop.preview); colormap('bone');
% drawnroi = images.roi.Circle('InteractionsAllowed', 'translate','Radius',R, 'LineWidth',1.5);
% draw(drawnroi);
% 
% 
% figure
% imagesc(VSDI.crop.preview); colormap('bone');
% drawncircle = drawcircle('Radius',6,'InteractionsAllowed', 'translate', 'LineWidth',1.5);


% drawnroi = drawcircle([],'Radius',40, 'InteractionsAllowed', 'translate', 'LineWidth',1.5);


%% DRAW CIRCULAR ROI WITH FUNCTION
clear
user_settings
nfish = 26;
VSDI = TORus('load',nfish);

VSDI.roi.labels_circ = {'dm4m_R', 'dm4m_L', ...
    'dm4_R','dm4_L',...
    'dm3_R','dm3_L',...
    'dm2_R','dm2_L',...
    'dldc_R','dldc_L',...
    'dldm_R','dldm_L',...
    'dm1_R','dm1_L'}

% for single hemisphere
VSDI.roi.labels_circ = {'dm4m_R', ...
    'dm4_R',...
    'dm3_R',...
    'dm2_R'};

VSDI.roi.labels_circ = {'dm4m_R', ... %FOR MAPS
    'dm4_R',...
    'dm2_R'};


close all
r = 4; % r = 4 for whole brain; 7 for 1x
[coord, mask] = roicir_draw(VSDI.crop.preview,VSDI.roi.labels_circ,r); 
[coord, mask] = roicir_draw(VSDI.backgr(:,:,VSDI.nonanidx(1)),VSDI.roi.labels_circ,r); 

circle.center = coord;
circle.R = r;
circle.mask = mask;

 roicirc_preview_multiple(VSDI.crop.preview, circle.center, circle.R); 
 roicirc_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(end)), circle.center, circle.R); 


imagesc(circle.mask(:,:,1))
axis image

VSDI.roi.circle.center = circle.center;
VSDI.roi.circle.R = circle.R;
VSDI.roi.circle.mask = circle.mask;
TORus('save',VSDI)

roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center, VSDI.roi.circle.R); 
% view rois from VSDI
 roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center, VSDI.roi.circle.R); 
ax1 = subplot(3,3,9);

%view one by one
for roi = 1:length(VSDI.roi.labels_circ)
     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(roi,:), VSDI.roi.circle.R); 
title([num2str(VSDI.ref) ':' VSDI.roi.labels_circ{roi}])
pause
end

%% ADD NEW CIRCULAR ROI TO THE ALREADY COLLECTION - !!! ONE BY ONE
clear
user_settings
nfish = 4;
VSDI = TORus('load',nfish);

newroi_label = {'dm2_L2'};

% newroi_label = {'dld_rR'};
newroi_row = numel(VSDI.roi.labels_circ)+1; %make sure it's an empty row

close all
r = VSDI.roi.circle.R;
[coord, mask] = roicir_draw(VSDI.crop.preview,newroi_label,r); 

circle.center = coord;
circle.R = r;
circle.mask = mask;

 roicirc_preview_multiple(VSDI.crop.preview, circle.center, circle.R); 
 roicirc_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(end)), circle.center, circle.R); 

VSDI.roi.labels_circ{newroi_row} =  newroi_label{1};
VSDI.roi.circle.center(newroi_row,:) = circle.center;
VSDI.roi.circle.mask(:,:,newroi_row) = circle.mask;
% TORus('save',VSDI)


% FINAL PREVIEW
roi_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center); 

%% DRAW RECTANGULAR CROP MASKS WITH FUNCTION
clear
user_settings
nfish = 23;
VSDI = TORus('load',nfish);

VSDI.roi.labels_rect = {'dm4_R'; 'dm2_R'};

close all
% [coord, mask] = roirect_draw(VSDI.crop.preview,VSDI.roi.labels_circ); 
[coord, mask] = roirect_draw(VSDI.backgr(:,:,VSDI.nonanidx(1)),VSDI.roi.labels_rect); 
[coord, mask] = roirect_draw(VSDI.backgr(:,:,VSDI.nonanidx(end)),VSDI.roi.labels_rect); 

roi_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(1)), coord); 
roi_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(end)), coord); 

imagesc(mask(:,:,1))
axis image

VSDI.roi.rect.coord = coord;
VSDI.roi.rect.mask = mask;
TORus('save',VSDI)

