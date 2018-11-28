function [fname_tetgen_femesh] = ...
    create_cells_femesh(fname_cell,fname_out_tetgen,include_box,box_gap,Rratio_nucleus,...
    cell_shape_name,hreq)

ndim = 3;

if (strcmp(cell_shape_name,'cylinders'))
    [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell] = ...
        create_cylinders_geometry(fname_cell,Rratio_nucleus);
elseif (strcmp(cell_shape_name,'ellipses'))
    [ncell,facets_cell,facets_labels_cell,nodes_cell,pt_in_cell] = ...
        create_ellipses_geometry(fname_cell,Rratio_nucleus);
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


%size(bd_attrib_allcell)
%disp('bd att');
if (include_box ~= 0)
    x1min = min(nodes(:,1)); x1max = max(nodes(:,1));
    x2min = min(nodes(:,2)); x2max = max(nodes(:,2));
    x3min = min(nodes(:,3)); x3max = max(nodes(:,3));
    
    x1len = x1max-x1min;
    x2len = x2max-x2min;
    x3len = x3max-x3min;
    
    x1gap = x1len*box_gap;
    x2gap = x2len*box_gap;
    x3gap = x3len*box_gap;
    
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
if (include_box ~= 0)
    regions = [regions;pt_in_box];
end
    
nregions = size(regions,1);

nodesindex = [1:nnodes]';

filename = [fname_out_tetgen,'.node'];
fid = fopen(filename,'w');
fprintf(fid, '%s\n', '# Part 1 - node list');
fprintf(fid, '%s\n', '# node count, 3 dim, no attribute, no boundary marker');
fprintf(fid, '%d %d %d %d\n', nnodes,ndim,0,0);
fprintf(fid, '%s\n', '# Node index, node coordinates');
if (fid ~= -1) 
  fprintf(fid, '%d %26.16f %26.16f %26.16f\n', [nodesindex,nodes]');
end
fclose(fid);

filename = [fname_out_tetgen,'.poly'];
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

tetgen_cmd = 'C:\Users\"Jing Rebecca LI"\WORK\WORK_RESEARCH\DMRI\CODE\gibbon\lib_ext\tetGen\win64\tetgen';
if (hreq > 0) 
    tetgen_options = ['-pqA','a',num2str(hreq)];
else
    tetgen_options = ['-pqA'];
end
disp(['Running command; ',tetgen_cmd,' ',tetgen_options,' ',filename]);
system([tetgen_cmd,' ',tetgen_options,' ',filename]);

[Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,...
          Nboundary,Ncmpt] = read_tetgen_new([fname_out_tetgen,'.1']);

fname_tetgen_femesh = [fname_out_tetgen,'.1'];
	
for icmpt = 1:Ncmpt
    Fac = [];
    for iboundary = 1:Nboundary
        Fac = [Fac,Fac_boundary_reorder{icmpt}{iboundary}];
    end
    [VOL_orig{icmpt}] ...
        = get_volume_mesh(Pts_cmpt_reorder{icmpt},Ele_cmpt_reorder{icmpt});
    [SA_orig{icmpt},SAu_orig{icmpt}] ...
        = get_surface_mesh_JRL(Pts_cmpt_reorder{icmpt},Fac);
end

VOL_total_orig = 0;
for icmpt = 1:Ncmpt
    VOL_total_orig  = VOL_total_orig + VOL_orig{icmpt};
end

for icmpt = 1:Ncmpt
    VOL_frac_orig{icmpt} = VOL_orig{icmpt}/VOL_total_orig;
end
for icmpt = 1:Ncmpt
    disp(['Tetgen: VF = ', num2str(VOL_frac_orig{icmpt}),', SA = ',num2str(SA_orig{icmpt}),...
        ', SAu = ',num2str(SAu_orig{icmpt})]);
end




