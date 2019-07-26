clear all;
format short

addpath SRC msh_files
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

fname_params_cells = 'params_cells_neuron_dendrites.in';
fname_params_simul_domain  = 'params_simul_domain_neuron_dendrites.in';
fname_params_simul_experi = 'params_simul_experi_neuron_dendrites_1direction.in';

DO_PLOTS = 1;

[params_cells,fname_cells] = create_geom(fname_params_cells);

[params_domain_geom,params_domain_pde,params_domain_femesh] ...
    = read_params_simul_domain(fname_params_simul_domain);

[DIFF_cmpts,kappa_bdys,IC_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,Ncmpt,Nboundary] ...
    = PREPARE_PDE(params_cells.ncell,params_cells.cell_shape,params_domain_geom,params_domain_pde);

fname_tetgen = [fname_cells];
if (params_cells.cell_shape == 3)
    fname_tetgen_femesh = fname_tetgen;
else
    [fname_tetgen_femesh] = ...
        create_femesh_fromcells(params_cells,fname_cells,params_domain_geom,params_domain_femesh,fname_tetgen);
end

ns = regexp(fname_tetgen_femesh,'/');
save_dir_name = fname_tetgen_femesh(ns(end)+1:end);
save_dir_path = ['saved_simul','/',save_dir_name];
tf = isfolder(save_dir_path);
if (~tf)
    call_cmd = ['mkdir ',save_dir_path];
    disp(call_cmd);
    eval([call_cmd]);
end

if (params_cells.cell_shape == 3)
    [mymesh] = read_mesh(fname_tetgen_femesh,DIFF_cmpts(1));
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

lscale = 4;
EigLim(OUT_cmpts_index) ...
    = DIFF_cmpts(OUT_cmpts_index).*(pi^2./(lscale.^2));
EigLim(IN_cmpts_index) ...
    = DIFF_cmpts(IN_cmpts_index).*(pi^2./(lscale.^2));
EigLim(ECS_cmpts_index) = DIFF_cmpts(ECS_cmpts_index).*(pi^2./(lscale.^2));

EigLim_min(OUT_cmpts_index) = -Inf;
EigLim_min(IN_cmpts_index) = -Inf;
EigLim_min(ECS_cmpts_index) = -Inf;

EigIntervals(1:Ncmpt-1) = 1;
EigIntervals(Ncmpt) = 30;

for icmpt = 1:Ncmpt
    fname = ['LAP_EIG_cmpt_',num2str(icmpt),...
        '_',num2str(EigLim_min(icmpt)),'_',num2str(EigLim(icmpt)/DIFF_cmpts(icmpt)),'.mat'];
    
    tf = isfile([save_dir_path,'/',fname]);
    if (tf)
        call_cmd = ['load ',save_dir_path,'/',fname];
        disp(call_cmd);
        eval([call_cmd]);
        
        EIG_value_cmpts{icmpt} = evalue;
        EIG_proj_cmpts{icmpt} = eproj;
        EIG_func_cmpts{icmpt} = func;
        ctime_laplace(1,icmpt) = ctime;
    else
        [EIG_value_cmpts{icmpt},EIG_proj_cmpts{icmpt},EIG_func_cmpts{icmpt},ctime_laplace(1,icmpt)]...
            = LAPLACE_EIG(mymesh.Pts_cmpt_reorder{icmpt},mymesh.Ele_cmpt_reorder{icmpt},...
            DIFF_cmpts(icmpt),VOL_cmpts(icmpt),EigLim(icmpt),EigLim_min(icmpt),EigIntervals(icmpt));
        
        LAP_EIG.evalue = EIG_value_cmpts{icmpt};
        LAP_EIG.eproj = EIG_proj_cmpts{icmpt};
        LAP_EIG.func = EIG_func_cmpts{icmpt};
        LAP_EIG.ctime = ctime_laplace(1,icmpt);
        
        call_cmd = ['save ',save_dir_path,'/',fname, ' ', '-struct', ' ', 'LAP_EIG'];
        disp(call_cmd);
        eval([call_cmd]);
    end
end

[EIG_length_cmpts] = MF_EIG_TO_LENGTH(EIG_value_cmpts,DIFF_cmpts);

experi_mf = experi_btpde;
nb = size(experi_mf.bvalues,2);
nexperi = size(experi_mf.bvalues,1);

[MF_JN_cmpts] = MF_JN(EIG_value_cmpts,DIFF_cmpts,experi_mf);

% Matrix Formalism effective diffusion tensor
ttt = cputime;

[DTENSOR_cmpts,DTENSOR_allcmpts] = MF_DTENSOR(VOL_cmpts,...
    IC_cmpts,EIG_value_cmpts,EIG_proj_cmpts,DIFF_cmpts,MF_JN_cmpts);

ctime_mf_dtensor = cputime - ttt;

if (DO_PLOTS ~= 0)
    cmpts_use = 1;
    PLOT_DTENSOR(DTENSOR_cmpts,DTENSOR_allcmpts,DIFF_cmpts);
    for nindex = 2
        diffdir = EIG_proj_cmpts{cmpts_use}(:,nindex);
        title_str = ['EigFun, l_s = ', num2str(EIG_length_cmpts{1}(nindex))...
                ', a_{1n} = [',num2str(diffdir'),']^T'];
        PLOT_PDESOLUTION(mymesh,EIG_func_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,...
            nindex,title_str);
    end
end


if (experi_common.ngdir_total == 1)
    if (~isempty(experi_btpde))
        % BTPDE
        [TOUT,YOUT,SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,difftime,ctime_btpde] ...
            = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts);
        [ADC_BTPDE_cmpts,ADC_BTPDE_allcmpts,ADC_allcmpts_S0] ...
            = FIT_SIGNAL(SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,experi_btpde.bvalues);
        [Sig_free,ADC_free_allcmpts] = ADCFREE(experi_btpde.bvalues,DIFF_cmpts,VOL_cmpts,IC_cmpts);
        

        % MF: Matrix Formalism
        [SIG_MF_cmpts,SIG_MF_allcmpts,ctime_mf] ...
            = SIG_MF(experi_mf,VOL_cmpts,IC_cmpts,EIG_value_cmpts,EIG_proj_cmpts);
        [ADC_MF_cmpts,ADC_MF_allcmpts,ADC_MF_allcmpts_S0] ...
            = FIT_SIGNAL(SIG_MF_cmpts,SIG_MF_allcmpts,experi_mf.bvalues);
        
        
        % MFGA: Matrix Formalism Gaussian Approximation
        [SIG_MFGA_cmpts,SIG_MFGA_allcmpts] ...
            = SIG_MFGA(experi_mf,VOL_cmpts,IC_cmpts,DTENSOR_cmpts);
        [ADC_MFGA_cmpts,ADC_MFGA_allcmpts,ADC_MFGA_allcmpts_S0] ...
            = FIT_SIGNAL(SIG_MFGA_cmpts,SIG_MFGA_allcmpts,experi_mf.bvalues);
        
        
        % PLOT BTPDE, MF and MFGA signals
        if (DO_PLOTS ~= 0)
            PLOT_SIGNAL(experi_btpde.bvalues,SIG_BTPDE_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_BTPDE_allcmpts,'BTPDE signal');
            PLOT_SIGNAL(experi_mf.bvalues,SIG_MF_allcmpts,Sig_free,ADC_MF_allcmpts_S0,ADC_MF_allcmpts,'MF signal');
            PLOT_SIGNAL(experi_mf.bvalues,SIG_MFGA_allcmpts,Sig_free,ADC_MFGA_allcmpts_S0,ADC_MFGA_allcmpts,'MFGA signal');
            PLOT_SIGNAL_BTPDE_MF(experi_btpde,SIG_BTPDE_allcmpts,SIG_MF_allcmpts,SIG_MFGA_allcmpts);
        end
    end
end