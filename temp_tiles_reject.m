    %% 6 - OVERLAID 12 TILES : from an average
    clear   
    user_settings
    
    set.manual_reject = 1;
    set.GSabsthres_reject =1;
    set.GSmethod_reject=1;
    set.force_include =1;
    
    
    pathsave = '/home/tamara/Documents/MATLAB/VSDI/TORus/plot/tiles/tiles_reject2';
    
   for  nfish = 8 %[2 3 4 8]
    
       
       VSDI = TORus('load', nfish);
%% SELECT CASES
cond_codes = unique(VSDI.condition(:,1));
cond_codes=  cond_codes(~isnan(cond_codes));

% cond_codes =[100 101 102 103 300 301 302 303];

%% SELECT EXCLUDED

rejectidx = [];

if set.manual_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
end

if set.GSabsthres_reject
    rejectidx = [rejectidx  makeRow(VSDI.reject.GSabs025)];

end

if set.GSmethod_reject
    rejectidx = [rejectidx makeRow(VSDI.reject.GSdeviat2sd)];
end 

if set.force_include
    rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
    
end

rejectidx = sort(unique(rejectidx));
    

%% 
    movie_ref = '_06filt3'; % input movie
    VSDmov = TORus('loadmovie',nfish,movie_ref);
    
    
    for ci = makeRow(cond_codes)

    sel_trials = find(VSDI.condition(:,1) == ci);
    sel_trials_clean = setdiff(sel_trials, rejectidx);
    
        
    %to plot average
    movie2plot = mean(VSDmov.data(:,:,:,sel_trials),4) ;
    movie2plot_clean = mean(VSDmov.data(:,:,:,sel_trials_clean),4) ;

    % settings
    tileset.start_ms = -12; % time in ms for first tile
    tileset.end_ms = 500%200;
    %           tileset.clims = [-0.9 0.9];
    tileset.clims = [-0.1 0.1];
    tileset.thresh = [-0.02 0.02];
    
    
    plot_tilemovie12frames(movie2plot, VSDI.timebase, tileset);
    titulo = [num2str(VSDI.ref) ,' AVERAGE cond', num2str(ci),'.plotthresh=',num2str(tileset.thresh(2)) ,'(all)'];
    sgtitle(titulo)
    
        name = [num2str(VSDI.ref),'.cond=',num2str(ci),'(all).jpg' ];
        saveas(gcf,fullfile(pathsave, name), 'jpg')
        close

    
    plot_tilemovie12frames(movie2plot_clean, VSDI.timebase, tileset);
    titulo2 = [num2str(VSDI.ref) ,' AVERAGE cond', num2str(ci),'.plotthresh=',num2str(tileset.thresh(2)) , '(clean)'];
    sgtitle(titulo2)
    
        name2 = [num2str(VSDI.ref), '.cond=',num2str(ci),'(clean).jpg' ];
        saveas(gcf,fullfile(pathsave, name2), 'jpg')
        close
        
%     BVmap= colormap_loadBV(); 
%     
%     mov = immovie(mat2im(movie2plot));
%         
%     implay(mov) 
    
    end %ci
    
   end %nfish 

   
