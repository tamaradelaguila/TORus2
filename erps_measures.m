%% PEAK-2-PEAK MEASSURE: 
% AMPLITUDE
% ONSET LATENCY


%% AMPLITUDE OF PEAK-2-PEAK

clear
user_settings
nfish = 8;%@ SET
VSDI = TORus('load',nfish);

refmovie = '_07filt4'; %@ SET
VSDroiTS = TORus('loadwave',nfish, refmovie);

timeser = VSDroiTS.filt407.data;

idx_trial = VSDI.nonanidx;
% idx_trial = find(VSDI.condition(:,1)==303);

idx_Stim = find_closest_timeidx(0, VSDI.timebase); 

for nroi = 1:length(VSDI.roi.labels)
for input.trialidx = 22 %find(VSDI.condition== cond_code);

cond_code=303;

input.data = squeeze(VSDroiTS.filt407.data(:,nroi,:));
input.timebase= VSDI.timebase; 
input.findwind_neg =[  -60;60 ];
input.findwind_pos = [  120;1000 ];
input.avew = 30; 
input.samplingtime = 6;
input.plot_on = 1;

[output, fig] = ana_erp_peak2peak(input);
title([num2str(VSDI.ref),':', VSDI.roi.labels{roi},'- condition=', num2str(cond_code)])
end %ntrial
end %nroi

%% FOR ERPs

for nroi = 1:length(VSDI.roi.labels)
        for cond_code= [300 301 303]

        input.trialidx = find(VSDI.condition== cond_code);

        input.data = squeeze(VSDroiTS.filt407.data(:,nroi,:));
        input.timebase= VSDI.timebase; 
        input.findwind_neg =[  -60;60 ];
        input.findwind_pos = [  120;1000 ];
        input.avew = 30; 
        input.plot_on = 1;

        [output, fig] = ana_erp_peak2peak(input);
        titulo = ['ERPs.', num2str(VSDI.ref),':', VSDI.roi.labels{nroi},'- condition=', num2str(cond_code)];
        title(titulo)
        plotpath = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/erps/measures';
        saveas(fig, fullfile(plotpath,titulo), 'jpg')
        close
        end
end %nroi