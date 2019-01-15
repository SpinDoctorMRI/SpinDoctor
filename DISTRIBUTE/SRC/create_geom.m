function [params_cells,fname_cells]=create_geom(fname_params_cells)

[params_cells,fname_cells] = read_params_cells(fname_params_cells);


if (params_cells.cell_shape == 1)
	create_ellipses_inputfile(params_cells.ncell,params_cells.Rmin,params_cells.Rmax,params_cells.dmin,params_cells.dmax,fname_cells);
elseif (params_cells.cell_shape == 2)
	create_cylinders_inputfile(params_cells.ncell,params_cells.Rmin,params_cells.Rmax,params_cells.dmin,params_cells.dmax,fname_cells,params_cells.Hcyl);
end