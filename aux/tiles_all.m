     %% 6 - OVERLAID 12 TILES : from an average
    clear
user_settings
nfish = 5; %@ SET

VSDI = TORus('load',nfish);

    movie_ref = '_06filt3'; % input movie
    VSDmov = TORus('loadmovie',nfish,movie_ref);
    
    
%      code2loop = {'100', '101', '102' , '103', '110', '111', '112', '113', '300', '301', '302','303' }; % @SET fish1
%      code2loop = {'200','201','202','203', '300', '301', '302','303' }; % @SET fish2
%      code2loop = {'100', '101', '102','103', '300', '301', '302','303' }; %SET fish3
%      code2loop = {'100', '101', '102','103'}; %SET fish4
     code2loop = {'1101', '2100', '2101','3100', '3101'}; %SET fish5
%      code2loop = {'1100', '1101','1401','2000', '2001', '2002', '2003', '3000', '3001', '3002', '4000', '4001','4002'}; %SET fish6
%      code2loop = {'1401'}; %SET fish6 - special case (control = 1100)

%%
     for ii =1:length(code2loop) %loop through all condition codes
        codeA = code2loop{ii};

        [idxA] = find(VSDI.condition(:,1)==str2num(codeA));
        %get control turning last number into 0
        temp = codeA; temp(end) = '0';
        [idxB] = find(VSDI.condition(:,1)==str2num(temp));

        %condition definition (adding mA only if there is one)
        mA = VSDI.list(idxA(1)).mA;
        if isnan(mA)
        cond_def = strcat('code#',codeA);
        else
        cond_def = strcat('code#',codeA,'(',num2str(VSDI.list(idxA(1)).mA),'mA)');
        end
    
     %to plot average
     movie2plot = mean(VSDmov.data(:,:,:,idxA),4, 'omitnan') ; 
    
     %average - control
     movie2plot_dif = mean(VSDmov.data(:,:,:,idxA),4, 'omitnan') - mean(VSDmov.data(:,:,:,idxB),4, 'omitnan') ; 
     movie2plot_dif(:,:,end) = movie2plot(:,:,end); %substitute the substracted background with a normal one
     
     % settings
          tileset.start_ms = -100; % time in ms for first tile
          tileset.end_ms = 800;
%           tileset.clims = [-0.9 0.9];
          tileset.clims =  [-0.6 0.6];
          tileset.thresh = [-0.2 0.2];
     
   % PLOT AND SAVE AVERAGE 

    plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
    titulo = [num2str(VSDI.ref) ,'AVE', cond_def];
    sgtitle(titulo)
   
                % Save and close

            path2save = fullfile('/home/tamara/Documents/MATLAB/VSDI/TORus//plot/tiles',num2str(VSDI.ref)) ; %@ SET
                if ~(isfolder(path2save))
                mkdir(path2save)
                end
            
            name2save = strcat(num2str(VSDI.ref),'AVE_',cond_def) ;  %@ SET (if you want)
               
            save_currentfig(path2save, name2save)
            close
            
    % PLOT AND SAVE AVERAGE MINUS CONTROL
    plot_tilemovie12frames(movie2plot_dif, VSDI.timebase, tileset);
    titulo2 = ['[',num2str(VSDI.ref) ,']','AVE minus CONTROL',cond_def];
    sgtitle(titulo2)
    
                    % Save and close
            name2save2 = strcat(num2str(VSDI.ref),'AVEdif_',cond_def) ;  %@ SET (if you want)
               
            save_currentfig(path2save, name2save2)
            close
            
     end