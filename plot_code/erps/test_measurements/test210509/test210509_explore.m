
%% GET MEASURES
% >>> RUN FIRST PART OF TEST210509


%% PLOT THRESHOLDS
climlocal = [min(noise_thresh.inf(:)), max(noise_thresh.sup(:))]
subplot(1,4,1); imagesc(noise_thresh.sup); set(gca, 'clim', climlocal)
subplot(1,4,2); imagesc(noise_thresh.inf); set(gca, 'clim', climlocal)
subplot(1,4,3); set(gca, 'clim', climlocal); colorbar
subplot(1,4,4); imagesc(noise_thresh.sup+noise_thresh.inf); set(gca, 'clim', climlocal)

colormap(jet)


climlocal = [min(noise_thresh.inf(:)), max(noise_thresh.sup(:))]
subplot(1,4,1); imagesc(noise_thresh.sup); set(gca, 'clim', climlocal)
subplot(1,4,2); imagesc(noise_thresh.inf); set(gca, 'clim', climlocal)
subplot(1,4,3); set(gca, 'clim', climlocal); colorbar
subplot(1,4,4); imagesc(noise_thresh.sup+noise_thresh.inf); set(gca, 'clim', climlocal)

%% PLOT PEAK-BASELINE RAW AND THRESHOLDED
%set to test
setbase =  [-500 0];
stdfactor = 1.5;
nplot = 2; 
        plot_peaklat = frames.peaklat;
        noise_thresh = return_noisethresh(movie, VSDI.timebase, setbase, stdfactor);
        
        outidx = frames.peakminusbasel < noise_thresh.sup & frames.peakminusbasel > noise_thresh.inf;
        
        plot_peakminusb = frames.peakminusbasel; 
        plot_peakminusb(outidx) = 0;
        
        climlocal = [min(plot_peakminusb(:)), max(plot_peakminusb(:))]*0.8;
figure
        subplot(1,3,1);  imagesc(frames.peakminusbasel(:,:,nplot)); set(gca, 'clim', climlocal); title('without thresh');
        subplot(1,3,2);  imagesc(plot_peakminusb(:,:,nplot)); set(gca, 'clim', climlocal); title('noise-thresh');
        subplot(1,3,3); set(gca, 'clim', climlocal); colorbar
        colormap(hot)
        
        sgtitle(['Peak-baseline. Noise thresh: baseline=',num2str(setbase(1)),'to',num2str(setbase(2)),'ms. std factor = ' num2str(stdfactor)])