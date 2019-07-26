function [fname_tetgen_femesh] = ...
   create_femesh_fromcells(params_cells,fname_cells,params_domain_geom,params_domain_femesh,fname_tetgen)

% create FE mesh from cells
% 
% Input:
%     1. params_cells is a structure with
%         a. 7 elements for spheres (cell_shape = 1):
%             cell_shape
%             ncell
%             Rmin
%             Rmax
%             dmin
%             dmax
%             para_deform 
%         b. 8 elements for cylinders (cell_shape = 2):
%             cell_shape
%             ncell
%             Rmin
%             Rmax
%             dmin
%             dmax
%             para_deform
%             Hcyl  
%
%     2. fname_cells
%
%     3. params_domain_geom is a structure with 3 elements:
%         Rratio_IN
%         include_ECS
%         ECS_gap
%
%     4. params_domain_femesh is a structure with 2 elements:
%         Htetgen
%         tetgen_cmd
%
%     5. fname_tetgen
% 
% Output: 
%     fname_tetgen_femesh

cell_shape = params_cells.cell_shape;
hreq = params_domain_femesh.Htetgen;
tetgen_cmd = params_domain_femesh.tetgen_cmd;
include_ECS = params_domain_geom.include_ECS;
ECS_gap = params_domain_geom.ECS_gap;
Rratio_IN = params_domain_geom.Rratio_IN;

ndim = 3;

if (cell_shape == 2)
    [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell,center,normal,Rcell,nslice_vec,nodes_ind_bottomring,nodes_ind_topring] ...
        = create_cylinders_geometry(fname_cells,Rratio_IN);    
elseif (cell_shape == 1)
    [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell,center,normal,Rcell] = ...
        create_ellipses_geometry(fname_cells,Rratio_IN);
elseif (cell_shape == 3)
    fname_tetgen_femesh = fname_cells;
    if (exist([fname_cells, '_elements.txt'], 'file') && exist([fname_cells, '_nodes.txt'], 'file')) == 0
        gmsh_to_fem(fname_cells);
    end
    return
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


nodes_allcell = nodes;
nnodes_allcell = size(nodes_allcell,1);

facets_allcell = facets;
nfacets_allcell = size(facets_allcell,2);

pt_in_allcell = pt_in_cell;
bd_attrib_allcell = boundary_attribute;

if (include_ECS == 2 & cell_shape == 2)
    theta = linspace(0,2*pi,21);
    thetavec = theta(1:end-1);
    xvec=cos(thetavec);
    yvec=sin(thetavec);
    
    nodes_all2d = [];
    for icell = 1:ncell                  
        center_one = center{icell};                   
        Rcell_one = Rcell{icell}*(1+ECS_gap);
        nodes_2d{icell}(:,1) = [xvec*Rcell_one(1,1)+center_one(1,1),center_one(1,1)]';
        nodes_2d{icell}(:,2) = [yvec*Rcell_one(1,1)+center_one(1,2),center_one(1,2)]';
        
        nodes_all2d = [nodes_all2d;nodes_2d{icell}];
        Rvec(icell) = Rcell_one(1,1);
    end
    
    Rmean = mean(Rvec(1:ncell));
    Rmin = min(Rvec(1:ncell));
    Rmax = max(Rvec(1:ncell));
    
    shp = alphaShape(nodes_all2d(:,1),nodes_all2d(:,2),Rmax*2);
    vv = area(shp);
    shp = alphaShape(nodes_all2d(:,1),nodes_all2d(:,2),Rmax,...
        'holeThreshold',vv/2,'regionThreshold',vv/2);
  
    if (numRegions(shp) > 1)
        disp('Can only have 1 region in ECS');
        stop
    end    
    [tri, xy] = shp.boundaryFacets();
    shp_bd = polyshape([xy(:,1)],[xy(:,2)]);
    clear shp_cell;
    for icell = 1:ncell
        ii = nodes_ind_topring{icell};
        nodes_topring{icell} = [nodes_cell{icell}(ii,1:3)];
         
        shp_cell(icell) = polyshape(nodes_topring{icell}(:,1),nodes_topring{icell}(:,2));

        ii = nodes_ind_bottomring{icell};
        nodes_bottomring{icell} = [nodes_cell{icell}(ii,1:3)];        

    end
    shp_ecs = subtract(shp_bd,union(shp_cell(1:ncell)));

    tri = shp_ecs.triangulation;
    % figure; 
    % triplot(tri.ConnectivityList,tri.Points(:,1),tri.Points(:,2));
    % title('ECS');
    % axis equal;
    
    Pts_ind = zeros(size(tri.Points,1),2);
    for icell = 1:ncell
        ii = nodes_ind_topring{icell};              
        for jj = 1:length(ii)
            kk = find(abs(nodes_cell{icell}(ii(jj),1)-tri.Points(:,1))<1e-6 & abs(nodes_cell{icell}(ii(jj),2)-tri.Points(:,2))<1e-6);  
            Pts_ind(kk,1:2) = [icell,ii(jj)];
        end      
    end
    ii = (find(Pts_ind(:,1)==0));
    Pts_ind(ii,1) = ncell+1;
    nodes_ecs = tri.Points(ii,1:2);
    Pts_ind(ii,2) = 1:length(ii);

    Pts_ecs_top = zeros(size(nodes_ecs,1),3);
    Pts_ecs_top(:,1:2) = [nodes_ecs];
    Pts_ecs_top(:,3) = nodes_cell{1}(nodes_ind_topring{1}(1),3);

    Pts_ecs_bottom = zeros(size(nodes_ecs,1),3);
    Pts_ecs_bottom(:,1:2) = [nodes_ecs];
    Pts_ecs_bottom(:,3) = nodes_cell{1}(nodes_ind_bottomring{1}(1),3);


    nodes = [nodes_allcell;Pts_ecs_bottom;Pts_ecs_top];

    Cmat = tri.ConnectivityList;
    Pts_top = zeros(size(tri.Points,1),3);
    Pts_top(:,1:2) = [tri.Points];
    Pts_top(:,3) = nodes_cell{1}(nodes_ind_topring{1}(1),3);

    Pts_bottom = zeros(size(tri.Points,1),3);
    Pts_bottom(:,1:2) = [tri.Points];
    Pts_bottom(:,3) = nodes_cell{1}(nodes_ind_bottomring{1}(1),3);

    Pts_both = [Pts_bottom;Pts_top];
    nfacets = size(Pts_bottom,1);
    Cmat_both = [Cmat;Cmat+nfacets];    

    
    tol = 1e-6;
    npt = size(Pts_both,1);
    rep_ind = zeros(npt,1);
    for ipt = 1:npt
        kk = find(abs(Pts_both(ipt,1)-nodes(:,1))<=tol & abs(Pts_both(ipt,2)...
                -nodes(:,2))<=tol & abs(Pts_both(ipt,3)-nodes(:,3))<=tol);
        rep_ind(ipt) = kk;
    end
    Cmat_new = zeros(size(Cmat_both));
    for ipt = 1:npt
        ii = find(Cmat_both==ipt);
        Cmat_new(ii) = rep_ind(ipt);        
    end    
    necs = size(nodes_ecs,1);
    [Cmat,Pts] = cylinder_connectivity(necs,necs);
    Cmat_out = Cmat+nnodes_allcell;

    
    nodes_box = [Pts_ecs_bottom; Pts_ecs_top];
    facets_box = [Cmat_new;Cmat_out']'-nnodes_allcell;

    pt_in_box = Pts_ecs_bottom(1,1:3);
    pt_in_box(1,3) = 0;

    bd_attrib_box = (nboundary+1)*ones(size(facets_box,2),1);

    nboundary = nboundary+1;
    
    nodes = [nodes_allcell;nodes_box];
    
    facets = [facets_allcell,facets_box+nnodes_allcell];

    bd_attrib = [bd_attrib_allcell; bd_attrib_box];    
  
elseif (include_ECS == 2 & cell_shape == 1)

    for icell = 1:ncell
        Rvec(icell) = mean(Rcell(icell,1:3));
    end
    
    Rmean = mean(Rvec(1:ncell));
    Rmin = min(Rvec(1:ncell));
    Rmax = max(Rvec(1:ncell));
        
    x1min = min(nodes(:,1)); x1max = max(nodes(:,1));
    x2min = min(nodes(:,2)); x2max = max(nodes(:,2));
    x3min = min(nodes(:,3)); x3max = max(nodes(:,3));  
    
    x1gap = Rmax*3;
    x2gap = Rmax*3;
    x3gap = Rmax*3;
    
    x1min = x1min-x1gap;
    x1max = x1max+x1gap;
    x2min = x2min-x2gap;
    x2max = x2max+x2gap;
    x3min = x3min-x3gap;
    x3max = x3max+x3gap;
    
    x1vec = linspace(x1min,x1max,80);
    x2vec = linspace(x2min,x2max,80);
    x3vec = linspace(x3min,x3max,80);
    
    [X1,X2,X3] = ndgrid(x1vec,x2vec,x3vec);
    
    xgap = ECS_gap*Rmax;
    
    distmat = Inf*ones(size(X1));
    for icell = 1:ncell
        dd = max(0,sqrt((X1-center{icell}(1,1)).^2+(X2-center{icell}(1,2)).^2+...
            +(X3-center{icell}(1,3)).^2)-Rvec(icell));
        distmat = min(distmat,dd);
    end
 
    [facets_box,nodes_box] = isosurface(X1,X2,X3,distmat,xgap);

    ii = find(distmat<=xgap*0.9 & distmat>=xgap*0.1);
    if (~isempty(ii))
        pt_in_box = [X1(ii(1)),X2(ii(1)),X3(ii(1))];
    else
        disp(['did not find pt in']);
        stop
    end
    
    facets_box = facets_box';

    %figure; h=trisurf(facets_box',nodes_box(:,1),nodes_box(:,2),nodes_box(:,3)); set(h,'facealpha',0.1); view(3); axis equal;
 
    bd_attrib_box = (nboundary+1)*ones(size(facets_box,2),1);
    nboundary = nboundary+1;
    
    nodes = [nodes_allcell;nodes_box];
    facets = [facets_allcell,facets_box+nnodes_allcell];
    bd_attrib = [bd_attrib_allcell; bd_attrib_box];

elseif (include_ECS == 1)
    x1min = min(nodes(:,1)); x1max = max(nodes(:,1));
    x2min = min(nodes(:,2)); x2max = max(nodes(:,2));
    x3min = min(nodes(:,3)); x3max = max(nodes(:,3));
    
    x1len = x1max-x1min;
    x2len = x2max-x2min;
    x3len = x3max-x3min;
    
    x1gap = x1len*ECS_gap;
    x2gap = x2len*ECS_gap;
    x3gap = x3len*ECS_gap;
    
    x1min = x1min-x1gap;
    x1max = x1max+x1gap;
    x2min = x2min-x2gap;
    x2max = x2max+x2gap;
    x3min = x3min-x3gap;
    x3max = x3max+x3gap;
    
    [nodes_box,facets_box,pt_in_box] = create_box_geometry(x1min,x1max,x2min,x2max,x3min,x3max);
    
    pt_in_box = [x1min+x1gap/2,x2min+x2gap/2,x3min+x3gap/2];
    
    bd_attrib_box = (nboundary+1)*ones(size(facets_box,2),1);
    nboundary = nboundary+1;
    
    nodes = [nodes_allcell;nodes_box];
    facets = [facets_allcell,facets_box+nnodes_allcell];
    bd_attrib = [bd_attrib_allcell; bd_attrib_box];
else
    nodes = [nodes_allcell];
    facets = [facets_allcell];
    bd_attrib = [bd_attrib_allcell];
end

nfacets = size(facets,2);
nnodes = size(nodes,1);

holes = [];
nholes = size(holes,1);

regions = [pt_in_allcell];
if (include_ECS ~= 0)
    regions = [regions;pt_in_box];
end
    
nregions = size(regions,1);

nodesindex = [1:nnodes]';

filename = [fname_tetgen,'.node'];
fid = fopen(filename,'w');
fprintf(fid, '%s\n', '# Part 1 - node list');
fprintf(fid, '%s\n', '# node count, 3 dim, no attribute, no boundary marker');
fprintf(fid, '%d %d %d %d\n', nnodes,ndim,0,0);
fprintf(fid, '%s\n', '# Node index, node coordinates');
if (fid ~= -1) 
  fprintf(fid, '%d %26.16f %26.16f %26.16f\n', [nodesindex,nodes]');
end
fclose(fid);

filename = [fname_tetgen,'.poly'];
fid = fopen(filename,'w');

fprintf(fid, '%s\n', '# Part 1 - node list');
fprintf(fid, '%s\n', '#  0 indicates the node list is stored in file .node');
fprintf(fid, '%d\n', 0);
fprintf(fid, '%s\n', '# Part 2 - facet list');
fprintf(fid, '%s\n', '# facet count, yes boundary marker');
fprintf(fid, '%d %d\n', nfacets,1);

fprintf(fid, '%s\n', '# Node index, node coordinates');
if (fid ~= -1) 
  for ii = 1:nfacets
    fprintf(fid, '%d %d %d\n', 1,0,bd_attrib(ii));  
    fprintf(fid, '%d ', 3);
    for jj = 1:3
      fprintf(fid, '%d ', facets(jj,ii));
    end
    fprintf(fid, '\n');
  end
end



%%%Part 3 - hole list
%%%Holes in the volume are specified by identifying a point inside each hole.
%%%  One line: <# of holes>
%%%  Following lines list # of holes:
%%%    <hole #> <x> <y> <z>
%%%    ...

fprintf(fid, '%s\n', ' # Part 3 - hole list');
fprintf(fid, '%d\n', nholes);
if (nholes>=1) 
	for ih = 1:nholes
		fprintf(fid, '%d %f %f %f\n', ih,holes(ih,1),holes(ih,2),holes(ih,3));
	end
end


fprintf(fid, '%s\n', ' # Part 4 - region list');
fprintf(fid, '%d\n', nregions);
if (nregions >= 1) 
	for ir = 1:nregions
		fprintf(fid, '%d %f %f %f %d %f\n', ir,regions(ir,1),regions(ir,2),regions(ir,3),ir,0.1);
	end
end
fclose(fid);

if (hreq > 0) 
    tetgen_options = ['-pqA','a',num2str(hreq)];
else
    tetgen_options = ['-pqA'];
end
disp(['Running command "',tetgen_cmd,' ',tetgen_options,' ',filename,'"']);
disp(['*****Start Tetgen ']);
system([tetgen_cmd,' ',tetgen_options,' ',filename]);
disp(['*****End Tetgen']);

fname_tetgen_femesh = [fname_tetgen,'.1'];