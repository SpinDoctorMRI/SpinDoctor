clear all; 
format short
addpath(genpath('SRC'));

fname_params_cells = 'params_cells.in';
fname_params_simul_domain  = 'params_simul_domain.in';
fname_params_simul_experi = 'params_simul_experi.in';

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
            ngdir_total = experi_btpde.ngdir_total;
            Ncmpt = mymesh.Ncmpt;
            nexperi = length(experi_btpde.sdeltavec);
            nb = size(experi_btpde.bvalues,2);
            
            % initialize matrices that contain diffusion signal
            SIG_BTPDE_cmpts_alldir = nan*ones(ngdir_total, Ncmpt, nexperi, nb);
            SIG_BTPDE_allcmpts_alldir = nan*ones(ngdir_total, nexperi, nb);
            ctime_alldir = nan*ones(ngdir_total, nexperi, nb);
            
            for ib = 1:nb
                % customized gradient directions. points_gdir: n x 3 matrix
                points_gdir = mygdirs(ngdir_total);
                
                % temporary experiment setting
                experi_tmp = experi_btpde;
                experi_tmp.qvalues = experi_btpde.qvalues(:,ib);
                experi_tmp.bvalues = experi_btpde.bvalues(:,ib);
                
                % compute signal
                [SIG_tmp_cmpts_alldir,SIG_tmp_allcmpts_alldir,ctime_tmp_alldir] ...
                            = SIG_BTPDE_HARDI(experi_tmp,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,points_gdir);
                
                SIG_BTPDE_cmpts_alldir(:, :, :, ib) = SIG_tmp_cmpts_alldir;
                SIG_BTPDE_allcmpts_alldir(:, :, ib) = SIG_tmp_allcmpts_alldir;
                ctime_alldir(:, :, ib) = ctime_tmp_alldir;
            end
            
            % plot signals
            if (0)
                for iexperi = 1:nexperi
                    for ib =1:nb
                        S0 = sum((IC_cmpts.*VOL_cmpts));
                        bv = experi_btpde.bvalues(iexperi,ib);
                        figure;
                        title_str = ['BTPDE SIGNAL. All Cmpts. ','Experi ', num2str(iexperi),', b= ',num2str(bv)];
                        plot3(points_gdir(:,1), points_gdir(:,2), squeeze(real(SIG_BTPDE_allcmpts_alldir(:,iexperi,ib)))/S0, 'k*');
                        title(title_str)
                        grid on
                    end
                end
            end
		end   
        %HADC
        if (~isempty(experi_hadc))
			[points,ADC_HADC_cmpts_alldir,ADC_HADC_allcmpts_alldir] ...
				= HADC_HARDI(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
            PLOT_HARDI(points,ADC_HADC_allcmpts_alldir);
        end
    end

end