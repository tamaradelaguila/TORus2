F =10;
[~, n] = colormap_loadBV2(9)
figure
s1= subplot(1,2,1);
imagesc(VSDI.backgr(:,:,10)); 
colormap(s1, colormap_loadBV2(n) )
colorbar

s2= subplot(1,2,2);
imagesc(VSDI.backgr(:,:,10)); 
colormap(s2, colormap_loadBV2(n*F) )
colorbar