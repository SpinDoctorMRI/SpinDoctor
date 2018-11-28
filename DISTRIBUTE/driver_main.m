clear; 

addpath SRC
addpath SRC/BTPDE SRC/DMRI SRC/FEM SRC/GEOMETRY SRC/TETGEN SRC/UTILITIES

SEQ_DEFINITIONS

fname_domain = 'InputFiles_Simulation\simulation_parameters_domain.in';
fname_experiment = 'InputFiles_Simulation\simulation_parameters_experiment.in';

[cell_shape,Rratio_nucleus,dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,ic_llimit,ic_ulimit,kappa_nc,kappa_ce,include_box,box_gap,...
    create_geom,fname_geom,ncell,Hcyl,Rmean,Rmin,Rmax,Htetgen] ...
    = read_simulation_parameters_domain(fname_domain);

[gdir,bvalues,qvalues,sdeltavec,bdeltavec,seqvec,npervec,rtol,atol,dt_out,snapshots,const_q,tetgen_cmd] ...
    = read_simulation_parameters_experiment(fname_experiment);

if (include_box == 1)
    box_str = 'box';
else
    box_str = 'nobox';
end
if (cell_shape == 1)
    cell_shape_name = 'ellipses';
elseif (cell_shape == 2)
    cell_shape_name = 'cylinders';
end
if (create_geom == 0)
    fname = fname_geom;
else
    fname = [cell_shape_name,num2str(ncell),'_R',num2str(Rmean)];
end
fname_cells_description = ['InputFiles_Geometry\',fname,'_description.in'];
if (create_geom == 0)
else
    if (cell_shape == 1)
        create_ellipses_inputfile(ncell,Rmean,Rmax,Rmin,fname_cells_description);
    else
        create_cylinders_inputfile(ncell,Rmean,Rmax,Rmin,Hcyl,fname_cells_description);
    end
end


fname_tetgen = ['InputFiles_Tetgen\',fname,'_',box_str];

[fname_tetgen_femesh] = create_cells_femesh(fname_cells_description,fname_tetgen,...
    include_box,box_gap,Rratio_nucleus,cell_shape_name,Htetgen,tetgen_cmd);
    
[TOUT,YOUT,MT,Ncmpt,Nboundary,Cell_cmpt,Box_cmpt,Nucleus_cmpt,nexperi,difftime,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,...
    DIFF_cmpts,IC_cmpts,UG] ...
    = solve_magnetization(fname_tetgen_femesh,ncell,Rratio_nucleus,...
    dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,ic_llimit,ic_ulimit,kappa_nc,kappa_ce,include_box,...
    gdir,qvalues,sdeltavec,bdeltavec,seqvec,npervec,rtol,atol);

[ADC_allcmpts,ADC_allcmpts_polydeg,ADC_allcmpts_S0,Deff_STA_allcmpts,...
    ADC,ADC_polydeg,ADC_S0,Deff_STA,MF_allcmpts,M0_allcmpts,S0_allcmpts,...
    MF,M0,S0,VOL,SA,SAu,VOL_frac,SoV]...
    = post_processing(MT,bvalues,Ncmpt,Nboundary,nexperi,sdeltavec,bdeltavec,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Fac_boundary_reorder,...
    DIFF_cmpts,UG);

plot_figures(UG,bvalues,difftime,...
    DIFF_cmpts,IC_cmpts,Ncmpt,ADC_allcmpts,ADC_allcmpts_S0,...
    Deff_STA_allcmpts,ADC,ADC_S0,...
    Deff_STA,MF_allcmpts,MF,VOL,SAu,VOL_frac);



    


    




