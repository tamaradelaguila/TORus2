nfish = 9;
VSDI = TORus('load', nfish ) ; 

n = length(VSDI.roi.labels_circ)
VSDI.roi.labels_circ{n+1,1} = 'dm4m_R2';
VSDI.roi.labels_circ{n+2,1} = 'dm2_R2';
VSDI.roi.labels_circ{n+3,1} = 'dm3_R2';
VSDI.roi.labels_circ{n+4,1} = 'dldm_R2';

VSDI.roi.circle.center(n+1,:) = [30 15]; %dm4
VSDI.roi.circle.center(n+2,:) = [38 36]; %dm2
VSDI.roi.circle.center(n+3,:) = [19 32]; %dm3
VSDI.roi.circle.center(n+4,:) = [10 40]; %dldm

% PREVISUALIZE AND ADJUST


                centers = VSDI.roi.circle.center(n+1:n+4, :) ;
                roicirc_preview_multiple(VSDI.crop.preview, centers, VSDI.roi.circle.R);
                
% ----------------------------------------------------------------
% MAKE NEW MASKS
% ----------------------------------------------------------------

xdim = size(VSDI.crop.preview,1); 
ydim = size(VSDI.crop.preview,2);


r= VSDI.roi.circle.R;
% Get the initial center position from the user
for roii = n+1:n+4
% Draw a cicle, and allow the user to change the position
coord = VSDI.roi.circle.center(roii,:)
drawncircle = drawcircle('Center', coord, 'InteractionsAllowed', 'translate', 'Radius', r );
VSDI.roi.circle.mask(:,:,roii) = createMask(drawncircle,xdim,ydim);

end

% CHECK CORRESPONDENCE BETWEEN ROIMASK - DRAWN ROI
roii = n+3
ax1 = subplot(1,2,1)
roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(roii,:), VSDI.roi.circle.R, ax1);

subplot(1,2,2)
imagesc(VSDI.roi.circle.mask(:,:,roii))
axis image

sgtitle(['roi:' VSDI.roi.labels_circ{roii}])



% ----------------------------------------------------------------
% SAVE VSDI WITH NEW DATA
% ----------------------------------------------------------------

% TORus('save', VSDI ) ; 

