function [Mrot] = rotate_nodes(u1,u2)
% Rotation matrix that rotates vector u1 to u2

u = u1+u2;
u = u/norm(u);
theta = pi;
ucross = cross(repmat(u,[1,3]),eye(3));
Mrot = cos(theta)*eye(3)+sin(theta)*ucross+(1-cos(theta))*u*u';

