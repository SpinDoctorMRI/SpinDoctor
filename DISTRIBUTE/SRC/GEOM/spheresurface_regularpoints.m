function [points,C,v]=spheresurface_regularpoints(r,N)

Ncount = 0;
a = 4*pi/N; d = sqrt(a);

Mnu = round(pi/d);

dnu = pi/Mnu;
dphi = a/dnu;

for m = 0:Mnu-1
    nu = pi*(m+0.5)/Mnu;
    Mphi = round(2*pi*sin(nu)/dphi);
    for n = 0:Mphi-1
        phi = 2*pi*n/Mphi;
        Ncount = Ncount +1;
        points(Ncount,1:3) = r*[sin(nu)*cos(phi),sin(nu)*sin(phi),cos(nu)];
    end
end
%figure; 
%plot3(points(:,1),points(:,2),points(:,3),':'); view(3); axis equal;

DT = delaunayTriangulation(points);
[C,v] = convexHull(DT);
%figure; 
%patch('faces',C,'vertices',points,'facealpha',0.1); view(3); axis equal;