%% s03 DEFINE ROIS
% Draw ROIs from brain areas.
clear
user_settings
nfish = 3;
VSDI = TORus('load',nfish);

% Draw ROIs and store in structure
VSDI.roi.labels = {'dm4_R', 'dm4_L',...
   'dm4m_R', 'dm4m_L',...
   'dm3_R', 'dm3_L',...
   'dm2_R', 'dm2_L',...
   'dldc_R', 'dldc_L',...
   'dldm_R', 'dldm_L',...
   'dm1_R', 'dmd1_L'}';

% DEFINE with the help of the function:
[manual_poly, manual_mask] = roi_draw(VSDI.crop.preview,VSDI.roi.labels); %in each column, one set of rois

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
TORus('save',VSDI);

roi_preview_multiple(VSDI.crop.preview, VSDI.roi.manual_poly); 

%% CREATE WAVES STRUCRURE (for extraction in '_extract_ROItimeseries):

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
nfish = 6;
VSDI = TORus('load',nfish);

close all
r = 4;
[coord, mask] = roicir_draw(VSDI.crop.preview,VSDI.roi.labels,r); 

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

% view rois from VSDI
 roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center, VSDI.roi.circle.R); 
ax1 = subplot(3,3,9)

%view one by one
for roi = 1:length(VSDI.roi.labels)
     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(roi,:), VSDI.roi.circle.R); 
title([num2str(VSDI.ref) ':' VSDI.roi.labels{roi}])
pause
end
