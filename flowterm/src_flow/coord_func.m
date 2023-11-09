function [x_g_pt_global,y_g_pt_global,z_g_pt_global] = coord_func(phi,coord)
%COORD_FUNC Summary of this function goes here
%   Detailed explanation goes here
x_g_pt_global = zeros(size(phi,4),size(phi,3));
y_g_pt_global = zeros(size(phi,4),size(phi,3));
z_g_pt_global = zeros(size(phi,4),size(phi,3));

for indp = 1:size(phi,3)
    for indele = 1:size(phi,4)
        x_g_pt_global(indele,indp) = sum(phi(1,:,indp,indele).*coord(1,:,indele));
        y_g_pt_global(indele,indp) = sum(phi(1,:,indp,indele).*coord(2,:,indele));
        z_g_pt_global(indele,indp) = sum(phi(1,:,indp,indele).*coord(3,:,indele));
    end
end

end

