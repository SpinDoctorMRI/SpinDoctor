clear all; 
format short

addpath SRC
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN
addpath Ex4

fname_params_cells = 'params_cells_Ex4.in';
fname_params_simul_domain  = 'params_simul_domain_Ex4.in';
fname_params_simul_experi = 'params_simul_experi_Ex4.in';

[params_cells,fname_cells] = create_geom(fname_params_cells);
	
PLOT_CELLS(params_cells.cell_shape,fname_cells);
	
[params_domain_geom,params_domain_pde,params_domain_femesh] ...
	= read_params_simul_domain(fname_params_simul_domain);

PLOT_SURFACE_TRIANGULATION(params_cells.cell_shape,fname_cells,params_domain_geom);

[DIFF_cmpts,kappa_bdys,IC_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,Ncmpt,Nboundary] ...
= PREPARE_PDE(params_cells.ncell,params_cells.cell_shape,params_domain_geom,params_domain_pde);

fname_tetgen = [fname_cells];
[fname_tetgen_femesh] = ...
   create_femesh_fromcells(params_cells,fname_cells,params_domain_geom,params_domain_femesh,fname_tetgen);
 
 
[mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);

if (~isempty(mymesh))
	PLOT_FEMESH(mymesh,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);

	[experi_common,experi_hadc,experi_btpde] ...
        = read_params_simul_experi(fname_params_simul_experi);
		
	[VOL_cmpts,SA_cmpts,SAu_cmpts,VOL_allcmpts,VF_cmpts,SoV_cmpts] ...
        = GET_VOL_SA(mymesh,experi_common.gdir);
		
	PLOT_GEOMETRY_INFO(cmpts_bdys_mat,OUT_cmpts_index,IN_cmpts_index,ECS_cmpts_index,VOL_cmpts,SA_cmpts,SAu_cmpts);
    
    if (experi_common.ngdir_total == 1) 
        % BTPDE
        if (~isempty(experi_btpde))
            [TOUT,YOUT,MF_cmpts,MF_allcmpts,difftime,BTPDE_elapsed_time] ...
                = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
            PLOT_MAGNETIZATION(mymesh,YOUT,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);
            [ADC_cmpts,ADC_allcmpts,ADC_allcmpts_S0] = FIT_SIGNAL(MF_cmpts,MF_allcmpts,experi_btpde.bvalues);       
            [Sig_free,ADC_free_allcmpts] = ADCFREE(experi_btpde.bvalues,DIFF_cmpts,VOL_cmpts,IC_cmpts);
            PLOT_SIGNAL(experi_btpde.bvalues,MF_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_allcmpts)
            PLOT_ADC(ADC_cmpts,ADC_allcmpts,DIFF_cmpts,'BTPDE');
            PLOT_TIMING(BTPDE_elapsed_time,mymesh);
        end
        %HADC
        if (~isempty(experi_hadc))
            [ADC_HADC_cmpts,ADC_HADC_allcmpts,HADC_elapsed_time] ...
                = HADC(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
            PLOT_ADC(ADC_HADC_cmpts,ADC_HADC_allcmpts,DIFF_cmpts,'HADC');
            PLOT_TIMING(HADC_elapsed_time,mymesh);
        end
        [ADC_STA_cmpts,ADC_STA_allcmpts] = STA(experi_common,DIFF_cmpts,VOL_cmpts,SAu_cmpts,IC_cmpts);	
		PLOT_ADC(ADC_STA_cmpts,ADC_STA_allcmpts,DIFF_cmpts,'STA');
        
    else
		% BTPDE
		if (~isempty(experi_btpde))	
			[SH_points,ADC_BT_cmpts_alldir,ADC_BT_allcmpts_alldir] ...
				= BTPDE_HARDI(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
            PLOT_HARDI(SH_points, ADC_BT_allcmpts_alldir);            
        end   
        %HADC
        if (~isempty(experi_hadc))
			[points,ADC_HADC_cmpts_alldir,ADC_HADC_allcmpts_alldir] ...
				= HADC_HARDI(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
            PLOT_HARDI(points,ADC_HADC_allcmpts_alldir);
        end
    end

end