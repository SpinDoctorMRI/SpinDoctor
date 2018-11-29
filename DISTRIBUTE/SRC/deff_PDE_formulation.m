function [deff_PDE_formulation_src,deff_PDE_formulation_src_time,ADC_PDE_formulation,Ncmpt,DIFF_cmpts] ...
    = deff_PDE_formulation(fname_tetgen_femesh,ncell,Rratio_nucleus,...
    dcoeff_nucleus,dcoeff_cytoplasm,dcoeff_exterior,...
    include_box,...
    gdir,sdeltavec,bdeltavec,seqvec,npervec,rtol,atol)

% diffusion equation (zero IC) to get the time-dependent diffusion coefficient

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG
global BDELTA SDELTA SEQ OGSEPER

disp(['In PDE formulation']);

tic

disp(['Reading tetgen Mesh ', fname_tetgen_femesh]);

[Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder,...
    Nboundary,Ncmpt] = read_tetgen_new([fname_tetgen_femesh]);

DIFF_cmpts = dcoeff_cytoplasm*ones(1,Ncmpt);
Cell_cmpt = 1:ncell;

if (include_box == 1)
    Box_cmpt = Ncmpt;
    Box_boundary = Nboundary;
    DIFF_cmpts(1,Box_cmpt) = dcoeff_exterior;
    
else
    Box_cmpt = [];
    Box_boundary = [];
end
if (Rratio_nucleus > 0)
    Nucleus_cmpt = ncell+1:2*ncell;
    Nucleus_boundary = ncell+1:2*ncell;
    DIFF_cmpts(1,Nucleus_cmpt) = dcoeff_nucleus;
 
else
    Nucleus_cmpt = [];
    Nucleus_boundary = [];    
end

UG = gdir';
UG = UG/norm(UG);

[model_FEM_matrices] = deff_PDE_formulation_FEMat(Ncmpt,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts);


toc

for icmpt = 1:Ncmpt
    [VOL(icmpt)] ...
        = get_volume_mesh(Pts_cmpt_reorder{icmpt},Ele_cmpt_reorder{icmpt});
end
for icmpt = 1:Ncmpt

    % The 6 FE matrices can be obtained in the following way.
    FEM_M = model_FEM_matrices{icmpt}.M;
    FEM_K = model_FEM_matrices{icmpt}.K;
    FEM_A = model_FEM_matrices{icmpt}.A;
    FEM_Q = model_FEM_matrices{icmpt}.Q;
    FEM_G = DIFF_cmpts(icmpt)*model_FEM_matrices{icmpt}.G;
     
   

    ODEsolve_atol = atol;
    ODEsolve_rtol = rtol;
    
    options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',...
        ODEsolve_rtol,'Vectorized','on','Stats','off',...
        'Jacobian',@odejac_bt_includeb);
    disp('DEFF PDE MODEL ***Uncoupled: start ode solver ode23t');
    nexperi = length(sdeltavec);
    for iexperi = 1:nexperi
        SDELTA = sdeltavec(iexperi);
        BDELTA = bdeltavec(iexperi);
        SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
        omega = 2*pi*npervec(iexperi)/SDELTA;
        OGSEPER = 1./omega*2*pi;%% set up number for OGSE
        QVAL = 0;
        TLIST = [0,SDELTA+BDELTA];
        ICC = zeros(size(FEM_M,1),1);
        sol = ode23t(@odefun_bt_includeb,TLIST,ICC,options);
        deff_PDE_formulation_src{iexperi}{icmpt} = FEM_G.'*sol.y/VOL(icmpt)/VOL(icmpt)/DIFF_cmpts(icmpt);
        deff_PDE_formulation_src_time{iexperi}{icmpt} = sol.x;         
        hvec = deff_PDE_formulation_src{iexperi}{icmpt};         
        tvec11 = deff_PDE_formulation_src_time{iexperi}{icmpt};
        Ftvec11 = seqintprofile(tvec11);
        a = trapz(tvec11,Ftvec11.*hvec*VOL(icmpt))/trapz(tvec11,Ftvec11.^2);
        ADC_PDE_formulation(icmpt,iexperi) = DIFF_cmpts(icmpt)*(1-a);
        
    end
end
