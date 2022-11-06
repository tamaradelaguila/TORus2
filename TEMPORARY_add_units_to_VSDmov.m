clear   
ref_movie = '_18filt6';
    
%     ref_movie = '_18filt6';
for nfish = [1:4 8:12]
   VSDmov = TORus('loadmovie', nfish, ref_movie);
   VSDmov.units = 'dF';
TORus('savemovie',VSDmov) 
end
blob()
% test wave 