function PLOT_HARDI_PT(points,SIG_allcmpts_alldir,fig_title)

% plot ngdir original directions in which the ADC or SIG was simulated
% 
% Input: 
%     1. points (ngdir directions)
%     2. ADC_allcmpts_alldir
%     3. fig_title
%
% Output:
%     1 figure with title of fig_title

ngdir = size(points,1);
DT = delaunayTriangulation(points);
[C_sph,v] = convexHull(DT);

SIG_interp = SIG_allcmpts_alldir;
SIG_alldir = SIG_interp.*points;
        
figure; hold on;
ADC_Ellipsoid = trisurf(C_sph,SIG_alldir(:,1),SIG_alldir(:,2),...
    SIG_alldir(:,3),SIG_interp,'facealpha',0.9); 

view([2,2,1]); 
caxis([0,1]); 
colorbar; 
axis equal;
grid on;
ADC_Ellipsoid.EdgeColor = 'none';

xlabel('x'); ylabel('y'); zlabel('z');
title(['SIG in ',num2str(ngdir),' directions. ',fig_title]);

plot3(SIG_allcmpts_alldir.*points(:,1),SIG_allcmpts_alldir.*points(:,2),...
    SIG_allcmpts_alldir.*points(:,3),'k.','markersize',5);