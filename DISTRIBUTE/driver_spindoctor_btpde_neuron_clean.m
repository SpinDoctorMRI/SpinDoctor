clear all;
format short

addpath SRC msh_files
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

fname_params_cells = 'params_cells_neuron_hardi.in';
fname_params_simul_domain  = 'params_simul_domain_neuron_hardi.in';
fname_params_simul_experi = 'params_simul_experi_neuron_hardi.in';
% fname_params_cells = 'params_cells_neuron.in';
% fname_params_simul_domain  = 'params_simul_domain_neuron.in';
% fname_params_simul_experi = 'params_simul_experi_neuron.in';

DO_PLOTS = 1;

% create cells if option 1 or 2
% 1 = ellipsoids, 2 = cylinders, 3 = neuron from msh file
[params_cells,fname_cells] = create_geom(fname_params_cells);

[params_domain_geom,params_domain_pde,params_domain_femesh] ...
    = read_params_simul_domain(fname_params_simul_domain);

[DIFF_cmpts,kappa_bdys,IC_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,Ncmpt,Nboundary] ...
    = PREPARE_PDE(params_cells.ncell,params_cells.cell_shape,params_domain_geom,params_domain_pde);

fname_tetgen = [fname_cells];

[fname_tetgen_femesh] = ...
    create_femesh_fromcells(params_cells,fname_cells,params_domain_geom,params_domain_femesh,fname_tetgen);

disp(['Reading FE mesh from:',fname_tetgen_femesh]);

if (params_cells.cell_shape == 3)
    params_cells.para_deform = [0,0]';
    [mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);
else
    [mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);
end

if (DO_PLOTS ~= 0)
    PLOT_FEMESH(mymesh,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);
end

[experi_common,experi_hadc,experi_btpde] ...
    = read_params_simul_experi(fname_params_simul_experi);

[VOL_cmpts,SA_cmpts,SAu_cmpts,VOL_allcmpts,VF_cmpts,SoV_cmpts] ...
    = GET_VOL_SA(mymesh,experi_common.gdir);

if (experi_common.ngdir_total == 1)
    if (~isempty(experi_btpde))
        % BTPDE
        [TOUT,YOUT,SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,difftime,ctime_btpde] ...
            = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
        [ADC_BTPDE_cmpts,ADC_BTPDE_allcmpts,ADC_allcmpts_S0] ...
            = FIT_SIGNAL(SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,experi_btpde.bvalues);
        [Sig_free,ADC_free_allcmpts] = ADCFREE(experi_btpde.bvalues,DIFF_cmpts,VOL_cmpts,IC_cmpts);
        
        if (DO_PLOTS ~= 0)
            PLOT_SIGNAL(experi_btpde.bvalues,SIG_BTPDE_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_BTPDE_allcmpts,'BTPDE');
        end
        
    end
else
    if (~isempty(experi_btpde))
        
        ngdir_total = experi_btpde.ngdir_total;
        [points_gdir,graddir_index,negii] = HARDI_PTS(ngdir_total);
        
        [SIG_BTPDE_cmpts_hardi,SIG_BTPDE_allcmpts_hardi,ctime_btpde_hardi] ...
            = SIG_BTPDE_HARDI(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,...
            points_gdir,graddir_index,negii);
        
        if (DO_PLOTS ~= 0)
            nb = size(experi_btpde.bvalues,2);
            nexperi = size(experi_btpde.bvalues,1);
            for iexperi = 1:nexperi
                for ib = 1:nb
                    S0 = sum((IC_cmpts.*VOL_cmpts));
                    bv = experi_btpde.bvalues(iexperi,ib);
                    title_str = ['BTPDE. All Cmpts. ','Experi ',...
                        num2str(iexperi),', b= ',num2str(bv)];
                    PLOT_HARDI_PT(points_gdir,...
                        squeeze(real(SIG_BTPDE_allcmpts_hardi(:,iexperi,ib)))/S0,title_str);
                end
            end
        end
    end
end

