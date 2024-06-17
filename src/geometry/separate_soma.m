function femesh_soma = separate_soma(femesh,swc,neighbours)
%%femesh_soma classifies the points of a femesh as being in the soma or dendrites.
% 
%   femesh: finite element mesh of the full cell.
%   swc: Swc class of the full cell
%   neighbours: array of tetrahedra adjacencies in the finite element mesh.
% 
%   femesh_soma: finite element mesh of soma


% Collect points from mesh and soma nodes from swc file
p = [femesh.points{:}];
npoint = size(p,2);
p = reshape(p,3,1,npoint);
soma_centres = swc.position_data(swc.type_data == 1,:)';
soma_radii = swc.radius_data(swc.type_data == 1)';

% Construct an inital crude mask for the soma elements as being within the
% soma balls
mask = squeeze(any(vecnorm(p - soma_centres,2,1) < 1.3*soma_radii,2));
e = [femesh.elements{:}];
p = [femesh.points{:}];
initial_dendrite_element_map = ~all(ismember(e,find(mask)),1);
initial_soma_elements = find(~initial_dendrite_element_map);

% Must make sure the final geometry is a valid mesh with a watertight
% surface. From the crude mask we take the connected component of the
% element closest to the origin.
g = get_femesh_network(neighbours);
[~, ~, centers_initial_soma] = get_volume_mesh(p, e(:,~initial_dendrite_element_map));
dist = vecnorm(centers_initial_soma,2,1);
[~,source_soma] = min(dist);
g = rmnode(g,find(initial_dendrite_element_map));
bins = conncomp(g);

% Final map of soma elements and points
soma_element_map = initial_soma_elements(bins == bins(source_soma));

% Once the mask is found, the mesh is extracted and information stored.
femesh_soma = initialise_soma(p,e,soma_element_map);

