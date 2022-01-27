
%% 


%% SETTINGS

condition.


%% LOAD DATA 
user_Settings

    ref_wave = 'circ_filt517';

    [VSDI] = TORus('load',nfish);
    
    if flagcirc
    roi1 = name2idx(roiname1, VSDI.roi.labels_circ);
    roi2 = name2idx(roiname2, VSDI.roi.labels_circ);
    else 
    roi1 = name2idx(roiname1, VSDI.roi.labels);
    roi2 = name2idx(roiname2, VSDI.roi.labels);

    end
    
        selroi = [roi1 roi2];

    
    VSDroiTS = TORus('loadwave',nfish);
    waves = VSDroiTS.(ref_wave).data; %@ SET
    

    
%% GET WAVE 
idx0 = find(VSDI.condition(:,1) == 404, 1,'first') ;
wave0 = squeeze(movies(30,25,1:end-1,idx0)); 

%% FFT ANALYSIS

Fs = 166;            % Sampling frequency                    
T = 1/Fs;             % Sampling period       
L = length(wave0);             % Length of signal
t = (0:L-1)*T;        % Time vector

% Compute the Fourier transform of the signal.
Y = fft(wave0);
% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f and plot the single-sided amplitude spectrum P1. The amplitudes are not exactly at 0.7 and 1, as expected, 
...because of the added noise. On average, longer signals produce better frequency approximations.

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

