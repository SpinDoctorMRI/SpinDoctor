function [femesh_dendrites,element_markers] = find_dendrites_femesh(femesh,femesh_soma,neighbours)
%%FIND_DENDRITES_FEMESH extracts the finite element meshes for the
%%dendrites.
% 
%   femesh: finite element mesh of the full cell
%   femesh_soma: finite element mesh of soma
%   neighbours: array of tetrahedra adjacencies in the finite element mesh.
% 
%   femesh_dendrites: cell array of finite element meshes for dendrites

g = get_femesh_network(neighbours);
g= rmnode(g,femesh_soma.element_map);
[bins,~] = conncomp(g);
N = size(neighbours,2);
mask = true([N,1]);
mask(femesh_soma.element_map) = false;
dendrite_elements = find(mask);

p = [femesh.points{:}];
e = [femesh.elements{:}];

[femesh_dendrites,element_markers] = initialise_dendrites(p,e,bins,dendrite_elements);