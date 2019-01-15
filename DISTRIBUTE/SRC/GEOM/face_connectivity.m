function [Cmat,center] = face_connectivity(nodes,ind_nodes,ind_center)

vec1 = ind_nodes;
vec2 = circshift(vec1,1);
vec3 = ones(size(vec1))*ind_center;
Cmat = [vec1;vec2;vec3];
center = zeros(1,3);
center(1,1) = mean(nodes(:,1));
center(1,2) = mean(nodes(:,2));
center(1,3) = mean(nodes(:,3));
