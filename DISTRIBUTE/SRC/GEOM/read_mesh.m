function [mymesh] = read_mesh(fname, sigma0)

% Read FE mesh from msh_files
%
% Input:
%     1. fname
%     2. sigma0  
%
% Output:
%     1. mymesh is a structure with 10 elements:
%         Nnode
%         Nele
%         Nface
%         Pts_cmpt_reorder
%         Ele_cmpt_reorder
%         Pts_ind
%         Pts_boundary_reorder
%         Fac_boundary_reorder
%         Nboundary
%         Ncmpt


% [Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder, Nboundary,Ncmpt]

disp(['Reading from neuron FE mesh from ', fname]);

mymodel = createpde();
elements = load([fname, '_elements.txt']);
nodes = load([fname, '_nodes.txt']);

geometryFromMesh(mymodel,nodes',elements');
% pdeplot3D(mymodel);

%%%%%%%%%%%%%
specifyCoefficients(mymodel,'m',0,'d',1,'c',sigma0,'a', 0,'f',0);
nface = mymodel.Geometry.NumFaces; %number of face
%     normal dot (c grad u) + q*u = g
for iface = 1:nface
    applyBoundaryCondition(mymodel,'neumann','face',iface,'g',1,'q',1,'Vectorized','on'); % 3-D geometry
end
%%%%%%%%%%%%%%

[elems2faces, faces2nodes] = get_faces(elements);
faces2elems=entryInWhichRows(elems2faces);
bindex=find(faces2elems(:, 2)==0);           %index of boundary facets - belong to one element only
boundary=faces2nodes(bindex, :);             %nodes of boundary facets
facets=sort(boundary, 2);


mymesh.Ncmpt = 1;
mymesh.Nboundary = 1;

mymesh.Nnode = size(nodes,1);
mymesh.Nele = size(elements,1);
mymesh.Nface = mymodel.Geometry.NumFaces;

mymesh.Pts_cmpt_reorder{1} = mymodel.Mesh.Nodes;
mymesh.Ele_cmpt_reorder{1} = mymodel.Mesh.Elements;
mymesh.Pts_ind{1} = [1:size(mymodel.Mesh.Nodes,2)]';

mymesh.Pts_boundary_reorder{1}{1} = [1:size(mymodel.Mesh.Nodes,2)]';
mymesh.Fac_boundary_reorder{1}{1} = facets';


