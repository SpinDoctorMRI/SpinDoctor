function [ADC,ADC_allcmpts,ADC_allcmpts_S0,TOUT,YOUT,MF_allcmpts,difftime,elapsed_time] ...
    = BTPDE(experiment,mm,DIFF_cmpts,kappa_vec,IC_cmpts)

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG 
global BDELTA SDELTA SEQ OGSEPER 
  
  bvalues = experiment.bvalues;
  qvalues = experiment.qvalues;
  
  gdir = experiment.gdir;
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;
rtol = experiment.rtol;
atol = experiment.atol;

yes = 1;  no = 0;
ODEsolve_atol = atol;
ODEsolve_rtol = rtol;


disp(['Magnetization: Solving Bloch-Torrey PDE']);

UG = gdir';
UG = UG/norm(UG);

if(mm.Ncmpt ==1 | abs(max(kappa_vec)) <= 1e-16)
    DO_COUPLING = no;
else
    DO_COUPLING = yes;
end


nexperi = length(sdeltavec);

%disp(['Simulating ',num2str(nexperi),' experiments']);

nb = size(qvalues,2);

elapsed_time=zeros(nb, nexperi);


    for icmpt = 1:mm.Ncmpt
        %disp(['Working on cmpt ', num2str(icmpt)]);
        coordinates = mm.Pts_cmpt_reorder{icmpt};
        elements = mm.Ele_cmpt_reorder{icmpt};
        facets = [];
        GX = -sqrt(-1)*UG*coordinates;
        FEM_MAT{icmpt}.Q = sparse(length(coordinates),length(coordinates));
        for iboundary = 1:mm.Nboundary
            neumann = mm.Fac_boundary_reorder{icmpt}{iboundary}';
            neumann_nodes=unique(neumann);
            coeffs_flux_matrix=zeros(max(neumann_nodes),1);
            coeffs_flux_matrix(neumann_nodes)=kappa_vec(iboundary);
            if ~isempty(neumann)
                FEM_MAT{icmpt}.Q = FEM_MAT{icmpt}.Q + flux_matrixP1_3D(neumann,coordinates',coeffs_flux_matrix);
            end;
        end;
        [FEM_MAT{icmpt}.K,volumes]=stiffness_matrixP1_3D(elements',coordinates',DIFF_cmpts(icmpt));
        FEM_MAT{icmpt}.M=mass_matrixP1_3D(elements',volumes);
        FEM_MAT{icmpt}.A=mass_matrixP1_3D(elements',volumes,GX');
    end

for icmpt=1:mm.Ncmpt
    % calculate the IC
    IC_Pts{icmpt} = IC_cmpts(icmpt)*ones(size(mm.Pts_cmpt_reorder{icmpt},2),1);
end

if (DO_COUPLING == yes)
    tic
    [FEMcouple_MAT,FEMcouple_ind0,FEMcouple_indf] ...
        = generate_FEM_coupling(FEM_MAT,mm.Ncmpt,mm.Nboundary,...
        mm.Pts_cmpt_reorder,mm.Ele_cmpt_reorder,mm.Pts_ind,mm.Pts_boundary_reorder,mm.Fac_boundary_reorder);
    %stop
    toc
    % IC for coupled case
    IC_couple=zeros(size(FEMcouple_MAT.M,1),1);
    for icmpt=1:mm.Ncmpt
        IC_couple(FEMcouple_ind0(icmpt):FEMcouple_indf(icmpt),1) = IC_Pts{icmpt};
    end
    
else
    FEMcouple_MAT = [];
    FEMcouple_ind0 = [];
    FEMcouple_indf = [];    
end

%% solve ODE
for iexperi = 1:nexperi
    
    SDELTA = sdeltavec(iexperi);
    BDELTA = bdeltavec(iexperi);
    TE = SDELTA+BDELTA;
    SEQ = seqvec(iexperi);% for choosing case PGSE, OGSEcos or OGSEsin
    omega = 2*pi*npervec(iexperi)/SDELTA;
    OGSEPER = 1./omega*2*pi;%% set up number for OGSE
    
    disp(['    Experiment ',num2str(iexperi)]);
    %disp([SDELTA,BDELTA,npervec(iexperi)]);
    
    TLIST = [0,SDELTA+BDELTA];
    for ib = 1:nb
        b_start_time = clock;
        % global variable setting QVAL for ODE time stepping
        QVAL = qvalues(iexperi,ib);        
        
		disp(['      qvalue ',num2str(QVAL,'%.1e')]);
        
        difftime(iexperi) = seqdifftime;
        
        %% Solving for case of coupling between compartments.
        
        if (DO_COUPLING == yes)
            
            FEM_M = FEMcouple_MAT.M;
            FEM_K = FEMcouple_MAT.K;
            FEM_A = FEMcouple_MAT.A; %*seqprofile(state.time)*QVAL;
            FEM_Q = FEMcouple_MAT.Q;
            
            FEM_G = sparse(zeros(size(FEM_M,1),1));
            
            options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                'Jacobian',@odejac_bt_includeb);            
            %disp('***Coupled: start ode solve ode23t'); tic
            sol= ode23t(@odefun_bt_includeb,TLIST,IC_couple,options);
            %disp('***Coupled: end ode solve ode23t'); toc
            for icmpt = 1:mm.Ncmpt
                YOUT{iexperi}{ib}{icmpt} = sol.y(FEMcouple_ind0(icmpt):FEMcouple_indf(icmpt),:);
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
            end
        else
            %% Solving for case of no coupling between compartments.            
            for icmpt = 1:mm.Ncmpt
                FEM_M = FEM_MAT{icmpt}.M;
                FEM_K = FEM_MAT{icmpt}.K;
                FEM_A = FEM_MAT{icmpt}.A;
                FEM_Q = FEM_MAT{icmpt}.Q;
                FEM_G = sparse(zeros(size(FEM_M,1),1));
                options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                    'Jacobian',@odejac_bt_includeb);
                %disp('***Uncoupled: start ode solver ode23t'); tic
                
                ICC = IC_Pts{icmpt};
                if (max(abs(ICC))<=1e-16)
                    sol.y = zeros(size(ICC,1),size(TLIST,2));
                    sol.x = TLIST;
                else
                    sol = ode23t(@odefun_bt_includeb,TLIST,ICC,options);
                end
                %disp('***Uncoupled: end ode solver ode23t'); toc
                YOUT{iexperi}{ib}{icmpt} = sol.y;
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
            end            
        end
        elapsed_time(ib, iexperi)=etime(clock, b_start_time);
    end
end

[ADC_allcmpts,ADC_allcmpts_polydeg,ADC_allcmpts_S0,...
    ADC,ADC_polydeg,ADC_S0,MF_allcmpts,M0_allcmpts,S0_allcmpts,MF,M0,S0]...
    = post_processing_magnetization(MT,bvalues);


