% EXPLORING SURFACE REPRESENTATION


ref_movie= '_17filt5';% '_17filt5' ; '_18filt6'

VSDmov = TORus('loadmovie', 11, ref_movie);

frame = VSDmov.data(:,:,end, 10); 

% TURN FRAME INTO SURFACE
for x = 1:size(frame,1)
    for y = 1:size(frame,2)
        xs(x,y) = x;
        ys(x,y) = y;
        zs(x,y) = frame(x,y);
    end
end

% PLOT SURFACE TO SEE
figure
s= surf(xs,ys,frame); colormap bone
%h= surf(xs,ys,frame, 'FaceAlpha',0.5)
s.EdgeColor = 'none';
view([-37.5 -30])

% INTERSECT WITH PLANE
[M,c] = contourf( xs,ys,frame, [5000 5000]);
c.LineWidth = 3;

%% TURN INTO COORDS  
[x,y,z] = C2xyz(M) ;
L = length(x);

% PLOT ON THE IMAGE
figure
    imagesc(frame); colormap bone
    hold on
    
for ii = 1:L
    xx= x{ii};
    yy = y{ii};
    plot(yy, xx, 'linewidth' , 1.5)
end


%%

% 1. MOVIE AVERAGE 

% 2. SET THRESHOLD-PLANE

% 3.SET  FRAMES(ms) TO INTERSECT AND GET FROM COLORMAP AS MANY COLOURS AS
% FRAMES

% 4 FOR EACH FRAME: 
% TURN INTO SURFACE
% - GET ALL CONTOURS



%% doesnt work
% % turn contour into polygon
% L = size(M,2); 
% 
% i = 1;
% c = 1;
% thresh = 2000; 
% while i<=L
%    n =  M(2,i);
%    Npoints(c)  = n;
%    Nidx(c) = i; 
%    i = i+n+1;
%    c = c+1;
% end

%%
% for ii = 1:length(Nidx)-1
%    
% curr= Nidx(ii)+1;
% post = Nidx(ii+1); post = post-1;
% coord = M(:,curr:post) ;    
% l = size(coord,2); 
% COORD{ii} = coord; 
% 
% end
% COORD{end+1} = M(:,Nidx(end)+1 : end);
% 
% clear COORD curr post coord ii

% plot(polyshape(COORD{1}'))
% plot(polyshape(COORD{2}'))
% plot(polyshape(COORD{3}'))
% 
% imagesc(frame); colormap bone
% hold on
% coord = COORD{1}';
% plot(coord(:,1), coord(:,2), 'color' , 'w', 'linewidth' , 1.5)
% plot(polyshape(coord))
%

