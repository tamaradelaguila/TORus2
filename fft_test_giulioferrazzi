VSDI = TORus('load', 10);
% VSDmov = TORus('loadmovie',10, '_03diffperc');

trials = find(VSDI.condition(:,1) == 2003); 
VSDmov = TORus('loadmovie',10, '_01registered');

A = VSDmov.data(:,:,:,trials);
A(isnan(A)) = 0;
image = mean(A,4);

IMAGE = fftshift(fft(image,[],3),3);

imagesc(log(squeeze(abs(IMAGE(30,:,:)))))

