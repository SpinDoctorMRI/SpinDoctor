function [points,theta,normals] = local_frustrum(n,r1,r2,h)
%LOCAL_FRUSTRUM Summary of this function goes here
%   Detailed explanation goes here
    [x,y] = fibonacci_lattice(n);
    theta = 2*pi.*x;
    points = zeros(3,n);normals = zeros(3,n);
    z  = h*y;
    r = (r1 + (r2 - r1).*z/h);
    points(1,:) = cos(theta).*r;
    points(2,:) = sin(theta).*r;
    points(3,:) = z;
    z_normal = (r1-r2)/h;
    normals(1,:) = cos(theta)./sqrt(1+(z_normal)^2);
    normals(2,:) = sin(theta)./sqrt(1+(z_normal)^2);
    normals(3,:) = normals(3,:) + z_normal;
end

