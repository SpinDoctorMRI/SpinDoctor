function [params_cells,fname_cells] = create_geom(fname_params_cells)

% create cells (canonical configuration)
% 
% Input:
%     fname_params_cells
%     
% Output: 
%     1. params_cells is a structure with
%         a. 7 elements for spheres (cell_shape = 1):
%             cell_shape
%             ncell
%             Rmin
%             Rmax
%             dmin
%             dmax
%             para_deform 
%
%         b. 8 elements for cylinders (cell_shape = 2):
%             cell_shape
%             ncell
%             Rmin
%             Rmax
%             dmin
%             dmax
%             para_deform
%             Hcyl   
%     2. fname_cells

[params_cells,fname_cells] = read_params_cells(fname_params_cells);

if (params_cells.cell_shape == 1)
	create_ellipses_inputfile(params_cells.ncell,params_cells.Rmin,params_cells.Rmax,params_cells.dmin,params_cells.dmax,fname_cells);
elseif (params_cells.cell_shape == 2)
	create_cylinders_inputfile(params_cells.ncell,params_cells.Rmin,params_cells.Rmax,params_cells.dmin,params_cells.dmax,fname_cells,params_cells.Hcyl);
end