function [cells_params,fname_cells]=create_cells_paramsetry(fname_cells_parameters)

[cells_params,fname_cells] = read_cells_parameters(fname_cells_parameters);


if (cells_params.cell_shape == 1)
	create_ellipses_inputfile(cells_params.ncell,cells_params.Rmin,cells_params.Rmax,cells_params.dmin,cells_params.dmax,fname_cells);
elseif (cells_params.cell_shape == 2)
	create_cylinders_inputfile(cells_params.ncell,cells_params.Rmin,cells_params.Rmax,cells_params.dmin,cells_params.dmax,fname_cells,cells_params.Hcyl);
end