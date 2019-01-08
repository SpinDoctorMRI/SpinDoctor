clear all; 
format short
addpath SRC
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOMETRY SRC/TETGEN SRC/UTILITIES SRC/COMPUTE_mesh_normals

fname_cells_parameters = 'cells_parameters.in';
fname_simulation_domain = 'simulation_domain_parameters.in';
fname_simulation_experiment = 'simulation_experiment_parameters.in';

[cells_params,fname_cells] = create_geometry(fname_cells_parameters);
	
PLOT_CELLS(cells_params.cell_shape,fname_cells);
	
[domain_geom,domain_pde,domain_femesh] ...
	= read_simulation_domain_parameters(fname_simulation_domain);

PLOT_SURFACE_TRIANGULATION(cells_params.cell_shape,fname_cells,domain_geom);

[DIFF_cmpts,kappa_vec,IC_cmpts,Cell_cmpt,Box_cmpt,Nucleus_cmpt,Ncmpt,Nboundary] ...
= prepare_PDE(cells_params.ncell,cells_params.cell_shape,domain_geom,domain_pde);

fname_tetgen = [fname_cells];
[fname_tetgen_femesh] = ...
   create_femesh_fromcells(cells_params,fname_cells,domain_geom,domain_femesh,fname_tetgen);
 
%fname_tetgen_femesh = 'current_geometry_15sph.1';
 
[mymesh,boundary_mat] = read_tetgen(fname_tetgen_femesh,cells_params.para_deform,Ncmpt,Nboundary);

if (~isempty(mymesh))


	
	PLOT_FEMESH(mymesh,Cell_cmpt,Box_cmpt,Nucleus_cmpt);

	[experiment_common,experiment_hadc,experiment_btpde] ...
		= read_simulation_experiment_parameters(fname_simulation_experiment);

		
	[VOL,SA,SAu,VOL_allcmpts,VOL_frac,SoV] ...
		= get_volume_surface(mymesh,experiment_common.gdir);
		
	PLOT_GEOMETRY_INFO(boundary_mat,Cell_cmpt,Nucleus_cmpt,Box_cmpt,VOL,SA,SAu);	


	if (~isempty(experiment_btpde))
		[ADC,ADC_allcmpts,ADC_allcmpts_S0,TOUT,YOUT,MF_allcmpts,difftime,BTPDE_elapsed_time] ...
			= BTPDE(experiment_btpde,mymesh,DIFF_cmpts,kappa_vec,IC_cmpts);
			
		PLOT_MAGNETIZATION(mymesh,YOUT,Cell_cmpt,Box_cmpt,Nucleus_cmpt)
		[Sig_free,ADC_free_allcmpts] = ADCFREE(experiment_btpde.bvalues,DIFF_cmpts,VOL,IC_cmpts);
		PLOT_SIGNAL(experiment_btpde.bvalues,MF_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_allcmpts)	
		
		PLOT_ADC(ADC,ADC_allcmpts,DIFF_cmpts,'BTPDE');
		%PLOT_TIMING(BTPDE_elapsed_time);
	end

	if (~isempty(experiment_hadc))
		[ADC_DE,ADC_DE_allcmpts,HADC_elapsed_time] ...
			= HADC(experiment_hadc,mymesh,DIFF_cmpts,IC_cmpts);
		PLOT_ADC(ADC_DE,ADC_DE_allcmpts,DIFF_cmpts,'HADC');	
		%PLOT_TIMING(HADC_elapsed_time);
	end


	[ADC_STA,ADC_STA_allcmpts] = STA(experiment_common,DIFF_cmpts,VOL,SAu,IC_cmpts);	
	PLOT_ADC(ADC_STA,ADC_STA_allcmpts,DIFF_cmpts,'STA');


end







