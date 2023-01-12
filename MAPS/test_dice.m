%% import coordinates for one fish and transform into image
clear all
close all
clf

nube1 = load('/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel_todos_dm4.csv');
tit = 'final-dm4-grupo';

area = nube1(:,2);
coor = nube1(:,3:4); 

scatter(coor(:,1), coor(:,2), 10, area(:,1), 'filled')

ymax = max(coor(:,1));
xmax = max(coor(:,2));

im = zeros(4, xmax,ymax);
imalpha =zeros(4, xmax,ymax);
for A = 1:4
    idx = find(area ==A);
    coorA = coor(idx,:);
    for ii = 1:length(coorA)
        xi = coorA(ii,1);
        yi = coorA(ii,2);
    imalpha(A,xi,yi) = 1;
        im(A,xi,yi) = 1;

    
%     im(A,xi,yi) = A;

    end
end

%% CONVERT LOGIC IMAGE INTO RGB
roimap(1,:) = [1 1 1]; %blue
roimap(2,:) = [0 .6 .8]; %blue
roimap(3,:) = [1 .6 0]; %orange
roimap(4,:) = [1 1 0];...yellow 
roimap(5,:) = [.6 .2 .95];...purple 


    imA = label2rgb (squeeze(imalpha(1,:,:)), roimap(2,:)); 
    imB = label2rgb (squeeze(imalpha(2,:,:)), roimap(3,:)); 
    imC = label2rgb (squeeze(imalpha(3,:,:)), roimap(4,:)); 
    imD = label2rgb (squeeze(imalpha(4,:,:)), roimap(5,:)); 
%% import coordinates for one fish and transform into image
clear all
close all
clf

nube1 = load('/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel221004_dm2.csv');
tit = 'final-221004-dm2';

area = nube1(:,2);
coor = nube1(:,3:4); 

scatter(coor(:,1), coor(:,2), 10, area(:,1), 'filled')

ymax = max(coor(:,1));
xmax = max(coor(:,2));

im = zeros(4, xmax,ymax);
imalpha =zeros(4, xmax,ymax);
for A = 1:4
    idx = find(area ==A);
    coorA = coor(idx,:);
    for ii = 1:length(coorA)
        xi = coorA(ii,1);
        yi = coorA(ii,2);
    imalpha(A,xi,yi) = 1;
        im(A,xi,yi) = 1;

    
%     im(A,xi,yi) = A;

    end
end

%% CONVERT LOGIC IMAGE INTO RGB
roimap(1,:) = [1 1 1]; %blue
roimap(2,:) = [0 .6 .8]; %blue
roimap(3,:) = [1 .6 0]; %orange
roimap(4,:) = [1 1 0];...yellow 
roimap(5,:) = [.6 .2 .95];...purple 


    imA = label2rgb (squeeze(imalpha(1,:,:)), roimap(2,:)); 
    imB = label2rgb (squeeze(imalpha(2,:,:)), roimap(3,:)); 
    imC = label2rgb (squeeze(imalpha(3,:,:)), roimap(4,:)); 
    imD = label2rgb (squeeze(imalpha(4,:,:)), roimap(5,:)); 

% %% PLOT OVERLAPPING MAKS - METHOD 1 FUSING INTO SINGLE IMAGE 
% % imagesc(squeeze(im(4,:,:))); 
% % 
% % imshowpair(im1,im2)
% % 
%         
%         
% imsum = double(imA)+double(imB)+double(imC)+double(imD);
% imsum = imsum/4;
% 
% imsum = rescale(imsum, 0, 225);
%  
% 
% % hsvImage = rgb2hsv(imsum);  %# Convert the image to HSV space
% % hsvImage(:,:,2) = 2;           %# Maximize the saturation
% % imsum = hsv2rgb(hsvImage);  %# Convert the image back to RGB space
% 
% imagesc(imsum)

%% PLOT OVERLAPPING MAKS - METHOD 2 OVERLAPPING IMAGES AXIS
 
% figure
clf
axH = subplot(1,2,1);
axH.Visible = 'off';

ax1 = axes;
imagesc(ax1, imA, 'alphadata',(mean(imA,3)>0)*1);
% ax1.Visible = 'off';
axis image
ax1.Visible = 'off';


ax2 = axes;
% imagesc(ax2,imAct,'alphadata',imAct>thresh);
% imalpha2 = imB>0;
imagesc(ax2,imB,'alphadata',(mean(imB,3)>0)*0.8);
ax2.Visible = 'off';
axis image


ax3 = axes;
% imagesc(ax2,imAct,'alphadata',imAct>thresh);
% imalpha2 = imB>0;
imagesc(ax3,imC,'alphadata',(mean(imC,3)>0)*0.4);
ax3.Visible = 'off';
axis image

ax4 = axes;
imagesc(ax4, imD, 'alphadata',(mean(imD,3)>0)*0.2);
% ax1.Visible = 'off';
axis image
ax4.Visible = 'off';

linkprop([axH ax1 ax2 ax3 ax4],'Position');


%% DICE INDEX FOR ONE SUBJECT
for A1 = 1:4
    for A2 = 1:4
    im1 = squeeze(im(A1,:,:));
    im2 = squeeze(im(A2,:,:));
    
    diceM(A1,A2) = dice(im1,im2);
%     inters =im1.*im2;  
%     dice2 = (2*sum(inters(:)))/(sum(im1(:))+sum(im2(:)));
clear im1 im2 
    end
end


%% PLOT DICE MATRIX
subplot(1,2,2)

imagesc(diceM)
colorbar
xticks([1:4])
xticklabels({'A','B','C','D'})
xlabel('positions')
yticks([1:4])
yticklabels({'A','B','C','D'})
ylabel('positions')
title('dice values')
axis image

sgtitle(tit)

diceM

%% NON-OVERLAPPING PLOT (for R)
dodgenube = nube1; 

% dodge points from each condition to improve visibility in scatterplot
dy = -0.1;
dx = -0.1;%% import coordinates for one fish and transform into image
clear all
close all
clf

nube1 = load('/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel221004_dm2.csv');
tit = 'final-221004-dm2';

area = nube1(:,2);
coor = nube1(:,3:4); 

scatter(coor(:,1), coor(:,2), 10, area(:,1), 'filled')

ymax = max(coor(:,1));
xmax = max(coor(:,2));

im = zeros(4, xmax,ymax);
imalpha =zeros(4, xmax,ymax);
for A = 1:4
    idx = find(area ==A);
    coorA = coor(idx,:);
    for ii = 1:length(coorA)
        xi = coorA(ii,1);
        yi = coorA(ii,2);
    imalpha(A,xi,yi) = 1;
        im(A,xi,yi) = 1;

    
%     im(A,xi,yi) = A;

    end
end

%% CONVERT LOGIC IMAGE INTO RGB
roimap(1,:) = [1 1 1]; %blue
roimap(2,:) = [0 .6 .8]; %blue
roimap(3,:) = [1 .6 0]; %orange
roimap(4,:) = [1 1 0];...yellow 
roimap(5,:) = [.6 .2 .95];...purple 


    imA = label2rgb (squeeze(imalpha(1,:,:)), roimap(2,:)); 
    imB = label2rgb (squeeze(imalpha(2,:,:)), roimap(3,:)); 
    imC = label2rgb (squeeze(imalpha(3,:,:)), roimap(4,:)); 
    imD = label2rgb (squeeze(imalpha(4,:,:)), roimap(5,:)); 

% %% PLOT OVERLAPPING MAKS - METHOD 1 FUSING INTO SINGLE IMAGE 
% % imagesc(squeeze(im(4,:,:))); 
% % 
% % imshowpair(im1,im2)
% % 
%         
%         
% imsum = double(imA)+double(imB)+double(imC)+double(imD);
% imsum = imsum/4;
% 
% imsum = rescale(imsum, 0, 225);
%  
% 
% % hsvImage = rgb2hsv(imsum);  %# Convert the image to HSV space
% % hsvImage(:,:,2) = 2;           %# Maximize the saturation
% % imsum = hsv2rgb(hsvImage);  %# Convert the image back to RGB space
% 
% imagesc(imsum)

%% PLOT OVERLAPPING MAKS - METHOD 2 OVERLAPPING IMAGES AXIS
 
% figure
clf
axH = subplot(1,2,1);
axH.Visible = 'off';

ax1 = axes;
imagesc(ax1, imA, 'alphadata',(mean(imA,3)>0)*1);
% ax1.Visible = 'off';
axis image
ax1.Visible = 'off';


ax2 = axes;
% imagesc(ax2,imAct,'alphadata',imAct>thresh);
% imalpha2 = imB>0;
imagesc(ax2,imB,'alphadata',(mean(imB,3)>0)*0.8);
ax2.Visible = 'off';
axis image


ax3 = axes;
% imagesc(ax2,imAct,'alphadata',imAct>thresh);
% imalpha2 = imB>0;
imagesc(ax3,imC,'alphadata',(mean(imC,3)>0)*0.4);
ax3.Visible = 'off';
axis image

ax4 = axes;
imagesc(ax4, imD, 'alphadata',(mean(imD,3)>0)*0.2);
% ax1.Visible = 'off';
axis image
ax4.Visible = 'off';

linkprop([axH ax1 ax2 ax3 ax4],'Position');


%% DICE INDEX FOR ONE SUBJECT
for A1 = 1:4
    for A2 = 1:4
    im1 = squeeze(im(A1,:,:));
    im2 = squeeze(im(A2,:,:));
    
    diceM(A1,A2) = dice(im1,im2);
%     inters =im1.*im2;  
%     dice2 = (2*sum(inters(:)))/(sum(im1(:))+sum(im2(:)));
clear im1 im2 
    end
end


%% PLOT DICE MATRIX
subplot(1,2,2)

imagesc(diceM)
colorbar
xticks([1:4])
xticklabels({'A','B','C','D'})
xlabel('positions')
yticks([1:4])
yticklabels({'A','B','C','D'})
ylabel('positions')
title('dice values')
axis image

sgtitle(tit)

diceM

%% NON-OVERLAPPING PLOT (for R)

clear all
close all

nube1 = load('/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel_todos_dm2.csv');
tit = 'final-todos-dm2';

dodgenube = nube1; 

step = 0.1;

% dodge points from each condition to improve visibility in scatterplot
dy = -step;
dx = -step;
idx = find(nube1(:,2)==1);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;
idx=[];

dy = +step;
dx = -step;
idx = find(nube1(:,2)==2);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;

dy = -step;
dx = +step;
idx = find(nube1(:,2)==3);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;

dy = +step;
dx = +step;
idx = find(nube1(:,2)==4);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;

outputfolder = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final';
outputname = ['dodge_', tit];
writematrix(dodgenube, fullfile(outputfolder, outputname))


%% PLOT OVERLAPPING MAKS - METHOD 2 OVERLAPPING IMAGES AXIS
 
% figure
clf
axH = subplot(1,2,1);
axH.Visible = 'off';

ax1 = axes;
imagesc(ax1, imA, 'alphadata',(mean(imA,3)>0)*1);
% ax1.Visible = 'off';
axis image
ax1.Visible = 'off';


ax2 = axes;
% imagesc(ax2,imAct,'alphadata',imAct>thresh);
% imalpha2 = imB>0;
imagesc(ax2,imB,'alphadata',(mean(imB,3)>0)*0.8);
ax2.Visible = 'off';
axis image


ax3 = axes;
% imagesc(ax2,imAct,'alphadata',imAct>thresh);
% imalpha2 = imB>0;
imagesc(ax3,imC,'alphadata',(mean(imC,3)>0)*0.4);
ax3.Visible = 'off';
axis image

ax4 = axes;
imagesc(ax4, imD, 'alphadata',(mean(imD,3)>0)*0.2);
% ax1.Visible = 'off';
axis image
ax4.Visible = 'off';

linkprop([axH ax1 ax2 ax3 ax4],'Position');


%% DICE INDEX FOR ONE SUBJECT
for A1 = 1:4
    for A2 = 1:4
    im1 = squeeze(im(A1,:,:));
    im2 = squeeze(im(A2,:,:));
    
    diceM(A1,A2) = dice(im1,im2);
%     inters =im1.*im2;  
%     dice2 = (2*sum(inters(:)))/(sum(im1(:))+sum(im2(:)));
clear im1 im2 
    end
end


%% PLOT DICE MATRIX
subplot(1,2,2)

imagesc(diceM)
colorbar
xticks([1:4])
xticklabels({'A','B','C','D'})
xlabel('positions')
yticks([1:4])
yticklabels({'A','B','C','D'})
ylabel('positions')
title('dice values')
axis image

sgtitle(tit)

diceM

%% NON-OVERLAPPING PLOT (for R)
dodgenube = nube1; 

% dodge points from each condition to improve visibility in scatterplot
dy = -0.1;
dx = -0.1;
idx = find(nube1(:,2)==1);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;
idx=[];

dy = +0.1;
dx = -0.1;
idx = find(nube1(:,2)==2);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;

dy = -0.1;
dx = +0.1;
idx = find(nube1(:,2)==3);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;

dy = +0.1;
dx = +0.1;
idx = find(nube1(:,2)==4);
dodgenube(idx,3)= nube1(idx,3)+dy;
dodgenube(idx,4)= nube1(idx,4)+dx;

outputfolder = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final';
outputname = ['dodge_', tit];
writematrix(dodgenube, fullfile(outputfolder, outputname))

