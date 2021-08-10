%% EXPORTA FRECUENCIA INSTANTÁTEA EN .CSV

% Este código
% (1) te exporta en csv y/o excel los tiempos de cada espiga y la frecuencia instantánea correspondiente a partir de un archivo '.mat' exportado desde spike. 
% (2) plotea la frecuencia instantánea de un intervalo
% (3) exporta en csv y/o excel la frecuencia instantánea para cada estímulo 

% Para ello, antes de exportarlo hay que crear los canales de eventos y nombrarlos de la siguiente manera: 
%   'spikes' - latidos
%   'stim' - llegada del estímulo

% Al exportar, marcar la opción 'Use source channel name in variable names' (pero desmarcando la anterior 'USe source name...'),
... y seleccionar en 'Layout options'> 'waveform and times'. 

% En el código sólo hay que rellenenar manualmente donde pone 'SETTINGS'. 

% SETTINGS: 
clear
carpeta.in = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike'; %poner la ruta de la carpeta que la contiene
carpeta.out = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike'; %donde quieres que se guarde el csv

nombre.input = '210325spike.mat';
nombre.output = 'test_save' ;

export.excel = 1; % '1' = sí; '0'=no
export.csv = 1; 
% FIN SETTINGS

% CARGAR EL ARCHIVO Y COPIAR EL CANAL 'spikes'
load(fullfile(carpeta.in, nombre.input))
clearvars -except spikes carpeta nombre export

% OBTENER LA MEDIDA
beats= spikes.times;
frec_inst = 1./(diff(beats));

% flipbeats= fliplr(beats);
% frec_inst = 1./(diff(flipbeats));
n = length(frec_inst); 

localoutput(:,1) = frec_inst; 
localoutput(:,2) = beats(2:end); 

% EXPORTAR
    %WRITE CSV
    if export.csv == 1
    filename = fullfile(carpeta.out,['spike_n_freq_',nombre.output,'.csv']);
    writematrix(localoutput,filename)
    disp(strcat(nombre.output, 'csv is ready in:', carpeta.out))

    end

    % WRITE EXCEL
    if export.excel == 1
    excelname = fullfile(carpeta.out,['spike_n_freq_',nombre.output,'.xlsx']);
    writematrix (localoutput,excelname)
    disp(strcat(nombre.output, [' excel file is ready in:' carpeta.out ]))
    end



%% SI QUIERES PLOTEAR UNA FRANJA TEMPORAL DEL ARCHIVO DE SPIKE

% SETTINGS
desde = 0; %en segundos
hasta = 7;

y_limites = [0.5 1.5]; %si quieres control sobre 
% FIN SETTINGS


beats = localoutput(:,2); frec = localoutput(:,1);

idx  = find(desde<beats & beats<hasta); 

plot(beats(idx),frec(idx), '.'); 

if  ~isempty('limites_y')
    ylim(y_limites);
end

%% PARA OBTENER LOS VALORES EN TORNO AL RO
% en tal caso hay que guardar otro canal con los 'eventos' de llegada del
% registro óptico, y llamarlo 'RO'
% SETTINGS: 
clear
carpeta.in = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike'; %poner la ruta de la carpeta que la contiene
carpeta.out = '/home/tamara/Documents/MATLAB/VSDI/TORus/data/dataspike'; %donde quieres que se guarde el csv

nombre.input = '210325spike.mat';
nombre.output = 'test_save' ;

w.pre= 3; % ventana pre-estímulo(en segundos) , de los que se va a tomar las espigas para hacer la media de frecuencia instantánea
w.post = 3;% ventana post-estímulo(en segundos) 
% w.out = 0; % ventana de seguridad a descartar (en segundos) si se quiere dejar el artefacto del estímulo fuera; si no hay artefacto, poner '0'

export.excel = 1; % '1' = sí; '0'=no
export.csv = 1; 
% FIN SETTINGS


% CARGAR EL ARCHIVO Y COPIAR LOS CANALES 'spikes' Y 'stim'
load(fullfile(carpeta.in, nombre.input))
beats= spikes.times;
stim= stim.times;
clearvars -except beats stim carpeta nombre export w

% descartar beats en la ventana de descarte de espigas
if size(stim, 1) > size (stim,2); stim = stim'; end

% CLEAN STIMULUS ARTIFACTS. If there is noise related to the stimulus, clean spikes in a time window after the stimulus before actually counting
beats_c = beats;
% for ii = 1:length(stim)
%     idx_discard= find(stim(ii)<beats & beats<stim(ii)+w.out );
%     beats_c(idx_discard) = NaN;
% end

% beats_c = beats_c(~isnan(beats_c));

frec_inst = 1./(diff(beats_c));
n = length(frec_inst); 

% CALCULATE PERCENTAGE OF BRADYCARDIC RESPONSE from inst freq
for ii= 1:length(stim)
    
pre = stim(ii) -w.pre;
post = stim(ii)+w.post; 

idx_pre = find(pre < beats & beats < stim(ii)); 
idx_post = find(stim(ii) < beats & beats < post); 

meanpre = mean(frec_inst(idx_pre));
meanpost = mean(frec_inst(idx_post));
perc = ((meanpre - meanpost) / meanpre )*100;

localoutput(ii,1) = stim(ii); %time of stimulus onset
localoutput(ii,2) = meanpre; %mean beats pre
localoutput(ii,3) = meanpost; %mean beats post
localoutput(ii,4) = perc; %percentage of brady from inst freq ( >0 means brady)

end

% EXPORTAR
    % WRITE CSV
    if export.csv == 1
    filename = fullfile(carpeta.out,['stim_percbrad_',nombre.output,'.csv']);
    writematrix(localoutput,filename)
    disp(strcat(nombre.output, 'csv is ready in:', carpeta.out))

    end

    % WRITE EXCEL
    if export.excel == 1
    excelname = fullfile(carpeta.out,['stim_percbrad_',nombre.output,'.xlsx']);
    writematrix (localoutput,excelname)
    disp(strcat(nombre.output, [' excel file is ready in:' carpeta.out ]))
    end
