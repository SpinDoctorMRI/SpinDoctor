function [femesh_soma,soma_element_map]= separate_soma(femesh,swc,neighbours,soma_mesh_path)
%%FEMESH_SOMA classifies the points of a femesh as being in the soma or dendrites.
% 
%   femesh: (struct) finite element mesh of the full cell.
%   swc: (Swc) swc object of the full cell
%   neighbours: [Nelement x 4] array of tetrahedra adjacencies in the finite element mesh.
%   soma_mesh_path: (string) path surface mesh of soma.
% 
%   femesh_soma: finite element mesh of soma


% Construct an inital crude mask for the soma vertices
if nargin == 4
   [f,v,~] = read_ply(soma_mesh_path);
    mask_points = [femesh.points{:}]';
    soma_center = mean(v)';
    disp("Separating soma")
    tic
    soma_mask = intriangulation(v,f,mask_points);
else
    p = [femesh.points{:}];
    npoint = size(p,2);
    p = reshape(p,3,1,npoint);
    soma_centres = swc.position_data(swc.type_data == 1,:)';
    soma_center = soma_centres(:,1);
    soma_radii = swc.radius_data(swc.type_data == 1)';
    soma_mask = squeeze(any(vecnorm(p - soma_centres,2,1) < 1.1*soma_radii,2));
    mask_points = squeeze(p)';
end

% Construct a separate mask for all vertices in the neurites/processes.
segment_mask = false(size(soma_mask));
eps = 0.001;
[~,segments] =  make_segment_list(swc.position_data,1.3*swc.radius_data,swc.connectivity_data,swc.type_data,1:size(swc.position_data,1) ,true);
for iseg = 1:length(segments)
    [in, on, ~, out_near] =intersect(segments(iseg),mask_points,eps);
    segment_mask = segment_mask | in | on |out_near;
end

% Include in the soma mask all points we cannot assign.
soma_mask = soma_mask | ~segment_mask;



% Find all elements with all four vertices belonging in the soma mask.
e = [femesh.elements{:}];
p = [femesh.points{:}];
initial_dendrite_element_map = ~all(ismember(e,find(soma_mask)),1);
initial_soma_elements = find(~initial_dendrite_element_map);

% From the crude mask we take the connected component of the element closest to the origin.
g = get_femesh_network(neighbours);
[~, ~, centers_initial_soma] = get_volume_mesh(p, e(:,~initial_dendrite_element_map));
size(centers_initial_soma)
size(soma_center)
dist = vecnorm(centers_initial_soma - soma_center,2,1);
[~,source_soma] = min(dist);
g = rmnode(g,find(initial_dendrite_element_map));
bins = conncomp(g);

% Final map of soma elements and points
soma_element_map = initial_soma_elements(bins == bins(source_soma));

% Once the mask is found, the mesh is extracted and information stored.
femesh_soma = initialise_soma(p,e,soma_element_map);

