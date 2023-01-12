%% import coordinates for one fish and transform into image
clearvars -except diceM
close all
clf
 
roi = 'dm2' ; % dm4 dm2 dm2_n3 (provisional)

switch roi
    case 'dm4'
        group{1} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel220706.csv';
        group{2} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel220712.csv';
        group{3} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel221007.csv';
        group{4} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel221004.csv';
    case 'dm2'
        group{1} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel220706_dm2.csv';
        group{2} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel220712_dm2.csv';
        group{3} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel221004_dm2.csv';
        group{4} = '/home/tamara/Documents/THESIS/05mapas_dice/nubes_final/nubepixel221007_dm2.csv';
end

for Sidx = 1:length(group)
nube1 = load(group{Sidx});
tit = ['group-' roi];

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

%     imA = label2rgb (squeeze(imalpha(1,:,:)), roimap(2,:)); 
%     imB = label2rgb (squeeze(imalpha(2,:,:)), roimap(3,:)); 
%     imC = label2rgb (squeeze(imalpha(3,:,:)), roimap(4,:)); 
%     imD = label2rgb (squeeze(imalpha(4,:,:)), roimap(5,:)); 

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


%% DICE INDEX FOR ALL SUBJECTS - MANUALLY RUN THE CODE ADJUSTING THE INTPUT CSV FILE AND THE SUBJECT IDX
for A1 = 1:4
    for A2 = 1:4
    im1 = squeeze(im(A1,:,:));
    im2 = squeeze(im(A2,:,:));
    
    diceM(A1,A2, Sidx) = dice(im1,im2);
%     inters =im1.*im2;  
%     dice2 = (2*sum(inters(:)))/(sum(im1(:))+sum(im2(:)));
clear im1 im2 
    end
end

end

% AVERAGE
diceM = squeeze(mean(diceM,3));

%% PLOT DICE MATRIX

imagesc(diceM)
colorbar
xticks([1:4])
xticklabels({'A','B','C','D'})
xlabel('positions')
yticks([1:4])
yticklabels({'A','B','C','D'})
ylabel('positions')
title(['dice values.' tit])
axis image


diceM