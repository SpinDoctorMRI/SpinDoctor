function femesh_soma =separate_soma_swc(femesh,swc,neighbours,soma_mesh_path)
[f,v,~] = read_ply(soma_mesh_path);
testpoints = [femesh.points{:}]';
[~,segments] =  make_segment_list(swc.position_data,1.3*swc.radius_data,swc.connectivity_data,swc.type_data,1:size(swc.position_data,1) ,true); % Changed true to false 20/02
soma_center = mean(v);
disp("Separating soma")
tic
soma_mask = intriangulation(v,f,testpoints);
segment_mask = false(size(soma_mask));
eps = 0.001;
for iseg = 1:length(segments)
[in, on, ~, out_near] =intersect(segments(iseg),testpoints,eps);
segment_mask = segment_mask | in | on |out_near;
end
soma_mask = soma_mask | ~segment_mask;

e = [femesh.elements{:}];
p = [femesh.points{:}];
initial_dendrite_element_map = ~all(ismember(e,find(soma_mask)),1);
initial_soma_elements = find(~initial_dendrite_element_map);
g = get_femesh_network(neighbours);
[~, ~, centers_initial_soma] = get_volume_mesh(p, e(:,~initial_dendrite_element_map));
dist = vecnorm(centers_initial_soma - soma_center',2,1);
[~,source_soma] = min(dist);
g = rmnode(g,find(initial_dendrite_element_map));
bins = conncomp(g);

% Final map of soma elements and points
soma_element_map = initial_soma_elements(bins == bins(source_soma));

% Once the mask is found, the mesh is extracted and information stored.
femesh_soma = initialise_soma(p,e,soma_element_map);