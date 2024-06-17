function g = get_femesh_network(neighbours)
% % GET_FEMESH_NETWORK returns a graph containing adjacency information of
% the tetrahedra of a finite element mesh
% 
%   neighbours: array of tetrahedra adjacencies in the finite element mesh.
% 
%   g: graph object with each node corresponding to an element of a finite
%       element mesh.


% Get index pairs of adjacent tetrahedra
index_pairs = cell(4,1);
for i = 2:5
indices = [neighbours(1,:);neighbours(i,:)]';
mask = all(indices > 0,2);
indices = indices(mask,:);
index_pairs{i-1} = indices;
end

% Create a sparse matrix to store information
N = size(index_pairs{1},1) + size(index_pairs{2},1) + size(index_pairs{3},1) + size(index_pairs{4},1);
i = [index_pairs{1}(:,1);index_pairs{2}(:,1);index_pairs{3}(:,1);index_pairs{4}(:,1)];
j = [index_pairs{1}(:,2);index_pairs{2}(:,2);index_pairs{3}(:,2);index_pairs{4}(:,2)];
v = ones(size(i));

adj = sparse(i,j,v,size(neighbours,2),size(neighbours,2));

% Output graph object
g = graph(adj);