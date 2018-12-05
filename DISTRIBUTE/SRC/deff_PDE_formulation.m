function [ADC_PDE_formulation] ...
    = deff_PDE_formulation(gdir,sdeltavec,bdeltavec,seqvec,npervec,rtol,atol,...
    Ncmpt,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts)

% diffusion equation (zero IC) to get the time-dependent diffusion coefficient

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG
global BDELTA SDELTA SEQ OGSEPER

disp(['In PDE formulation']);


UG = gdir';
UG = UG/norm(UG);

[model_FEM_matrices] = deff_PDE_formulation_FEMat(Ncmpt,Pts_cmpt_reorder,Ele_cmpt_reorder,DIFF_cmpts);

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
       
        %deff_PDE_formulation_src{iexperi}{icmpt} = FEM_G.'*sol.y/VOL(icmpt)/VOL(icmpt)/DIFF_cmpts(icmpt);
        deff_PDE_formulation_src{iexperi}{icmpt} = FEM_G.'*sol.y/VOL(icmpt)/VOL(icmpt);

        deff_PDE_formulation_src_time{iexperi}{icmpt} = sol.x;

        hvec = deff_PDE_formulation_src{iexperi}{icmpt};
        tvec11 = deff_PDE_formulation_src_time{iexperi}{icmpt};
        Ftvec11 = seqintprofile(tvec11);
        a = trapz(tvec11,Ftvec11.*hvec*VOL(icmpt))/trapz(tvec11,Ftvec11.^2);
        
        ADC_PDE_formulation(icmpt,iexperi) = DIFF_cmpts(icmpt)-a;

        %ADC_PDE_formulation(icmpt,iexperi) = DIFF_cmpts(icmpt)*(1-a);
        %deff_PDE_formulation_src{iexperi}{icmpt} = FEM_G.'*sol.y/VOL(icmpt)/VOL(icmpt)/DIFF_cmpts(icmpt);
    end
end
