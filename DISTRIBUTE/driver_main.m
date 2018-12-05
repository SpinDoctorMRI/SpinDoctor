clear; 

addpath SRC
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOMETRY SRC/TETGEN SRC/UTILITIES

mydefinitions;
ncolor = length(colorvec_cell);

SEQ_DEFINITIONS

fname_domain = 'InputFiles_Simulation\simulation_parameters_domain.in';
fname_experiment = 'InputFiles_Simulation\simulation_parameters_experiment.in';

[cell_shape,Rratio_nucleus,dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,ic_llimit,ic_ulimit,kappa_nc,kappa_ce,include_box,box_gap,...
    create_geom,fname_geom,ncell,Hcyl,Rmean,Rmin,Rmax,Htetgen,para_deform] ...
    = read_simulation_parameters_domain(fname_domain);

[gdir,bvalues,qvalues,sdeltavec,bdeltavec,seqvec,npervec,...
    rtol_bt,atol_bt,rtol_deff,atol_deff,const_q,tetgen_cmd] ...
    = read_simulation_parameters_experiment(fname_experiment);

if (include_box ~= 0)
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

disp(['Reading tetgen Mesh ', fname_tetgen_femesh]);

[Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,...
    Nboundary,Ncmpt] = read_tetgen_new(fname_tetgen_femesh,[0,0]);

[VOL_orig,SA_orig,SAu_orig,VOL_allcmpts_orig,VOL_frac_orig,SoV_orig] = get_volume_surface(Nboundary,Ncmpt,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Fac_boundary_reorder,eye(3));

for icmpt = 1:Ncmpt
    disp(['No Deform Mesh: VF = ', mynum2str(VOL_frac_orig(icmpt)),...
        ', SoV = ', mynum2str(SoV_orig(icmpt)),', SA = ',mynum2str(SA_orig(icmpt)),...
        ', SAu = ',mynum2str(SAu_orig(icmpt,:))]);
end

disp(['Reading tetgen Mesh ', fname_tetgen_femesh]);
[Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,...
    Nboundary,Ncmpt] = read_tetgen_new([fname_tetgen_femesh],para_deform);

boundary_mat = zeros(Ncmpt,Nboundary);
for icmpt = 1:Ncmpt
    for ibd = 1:Nboundary
        boundary_mat(icmpt,ibd)=~isempty(Pts_boundary_reorder{icmpt}{ibd});
    end
end

% [VOL,SA,SAu_test,VOL_allcmpts,VOL_frac,SoV_test] ...
%     = get_volume_surface(Nboundary,Ncmpt,...
%     Pts_cmpt_reorder,Ele_cmpt_reorder,Fac_boundary_reorder,eye(3));
% UG = gdir';
% UG = UG/norm(UG);
% SAu2 = SAu_test*(UG'.^2);
% SoV2 = SoV_test*(UG'.^2);
% for icmpt = 1:Ncmpt
%     disp(['Deform Mesh: VF = ', mynum2str(VOL_frac(icmpt)),...
%         ', SoV2 = ', mynum2str(SoV2(icmpt)),', SA = ',mynum2str(SA(icmpt)),...
%         ', SAu2 = ',mynum2str(SAu2(icmpt,:))]);
% end

[VOL,SA,SAu,VOL_allcmpts,VOL_frac,SoV] ...
    = get_volume_surface(Nboundary,Ncmpt,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Fac_boundary_reorder,gdir);

for icmpt = 1:Ncmpt
    disp(['Deform Mesh: VF = ', mynum2str(VOL_frac(icmpt)),...
        ', SoV = ', mynum2str(SoV(icmpt)),', SA = ',mynum2str(SA(icmpt)),...
        ', SAu = ',mynum2str(SAu(icmpt,:))]);
end

[DIFF_cmpts,kappa_vec,IC_cmpts,Cell_cmpt,Box_cmpt,Nucleus_cmpt] ...
    = assign_diffkappa(Ncmpt,Nboundary,ncell,Rratio_nucleus,...
    dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    ic_nucleus,ic_cytoplasm,ic_exterior,kappa_nc,kappa_ce,...
    include_box,cell_shape_name);

[ADC_PDE_formulation] ...
    = deff_PDE_formulation(gdir,sdeltavec,bdeltavec,seqvec,npervec,rtol_deff,atol_deff,...
    Ncmpt,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts);

 [TOUT,YOUT,MT,difftime] ...
    = solve_magnetization(gdir,qvalues,sdeltavec,bdeltavec,seqvec,npervec,rtol_bt,atol_bt,...
    Ncmpt,Nboundary,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts,...
    Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,kappa_vec,IC_cmpts,ic_llimit,ic_ulimit);

[ADC_allcmpts,ADC_allcmpts_polydeg,ADC_allcmpts_S0,...
    ADC,ADC_polydeg,ADC_S0,MF_allcmpts,M0_allcmpts,S0_allcmpts,...
    MF,M0,S0] = post_processing(MT,bvalues);

UG = gdir';
UG = UG/norm(UG);

nexperi = length(MT);
nb = length(MT{1});

ADC_PDE_allcmpts = nan*ones(nexperi,1);
for iexperi = 1:nexperi
    ADC_PDE_allcmpts(iexperi,1) = sum((IC_cmpts.*VOL_frac)'.*ADC_PDE_formulation(:,iexperi))./sum((IC_cmpts.*VOL_frac)');
end

Sig_free = zeros(size(bvalues(:)));
vol = 0;
for icmpt = 1:Ncmpt
    Sig_free = Sig_free+IC_cmpts(1,icmpt)*VOL_frac(icmpt)*exp(-DIFF_cmpts(icmpt)*bvalues(:));
    vol = vol+IC_cmpts(1,icmpt)*VOL_frac(icmpt);
end
Sig_free = Sig_free/vol;
ADC_free_allcmpts = sum((IC_cmpts.*VOL_frac)'.*DIFF_cmpts')./sum((IC_cmpts.*VOL_frac)');

%% Deff short time approximation
ADC_STA = zeros(Ncmpt,nexperi);
for iexperi = 1:nexperi
    for icmpt = 1:Ncmpt
        [ADC_STA(icmpt,iexperi)] = deff_sta(DIFF_cmpts(icmpt),...
            VOL(icmpt),SAu(icmpt),sdeltavec(iexperi),bdeltavec(iexperi),...
            seqvec(iexperi),npervec(iexperi)); % short time approximation        
    end
end
ADC_STA_allcmpts = nan*ones(nexperi,1);
for iexperi = 1:nexperi
    ADC_STA_allcmpts(iexperi,1) = sum((IC_cmpts.*VOL_frac)'.*ADC_STA(:,iexperi))./sum((IC_cmpts.*VOL_frac)');
end

plot_results;
