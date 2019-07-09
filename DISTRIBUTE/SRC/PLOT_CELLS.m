function PLOT_CELLS(cell_shape,fname_cells)

% plot cells of canonical configuration
% 
% Input:
%     1. cell_shape
%     2. fname_cells
%
% Output: 
%     1 figure with title of "Cells of canonical configuration"

if (cell_shape == 2)
    [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell,center,normal,Rcell,nslice_vec,nodes_ind_bottomring,nodes_ind_topring] ...
        = create_cylinders_geometry(fname_cells,0);    
elseif (cell_shape == 1)
    [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell,center,normal,Rcell] = ...
        create_ellipses_geometry(fname_cells,0);
end

facets = [];
nodes = [];
offset = 0;

boundary_attribute = [];
nboundary = 0;
for icell = 1:ncell
	facets = [facets,facets_cell{icell}+offset];	
	nfacets = size(facets_cell{icell},2);
	nnodes = size(nodes_cell{icell},1);	        
	boundary_attribute = [boundary_attribute;facets_labels_cell{icell}'+nboundary];	
    nboundary = nboundary + max(facets_labels_cell{icell}');
	nodes = [nodes;nodes_cell{icell}];
	offset = offset+nnodes;
end

figure; 
h = trisurf(facets',nodes(:,1),nodes(:,2),nodes(:,3)); 
set(h,'facealpha',1); view(3); axis equal;
title('Cells of canonical configuration');