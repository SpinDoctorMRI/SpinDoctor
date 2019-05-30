function [TOUT,YOUT,MF_cmpts,MF_allcmpts,difftime,elapsed_time] ...
    = BTPDE(experiment,mymesh,DIFF_cmpts,T2_cmpts,kappa_bdys,IC_cmpts)

% solve Bloch-Torrey equation
% 
% Input:
%     1. experiment is a structure with 10 elements:
%         ngdir_total 
%         gdir        
%         sdeltavec   
%         bdeltavec   
%         seqvec       
%         npervec     
%         rtol        
%         atol        
%         qvalues     
%         bvalues        
%     2. mymesh is a structure with 10 elements:
%         Nnode
%         Nele
%         Nface
%         Pts_cmpt_reorder
%         Ele_cmpt_reorder
%         Pts_ind
%         Pts_boundary_reorder
%         Fac_boundary_reorder
%         Nboundary
%         Ncmpt
%     3. DIFF_cmpts
%     4. kappa_bdys
%     5. IC_cmpts
%
% Output:
%     1. TOUT
%     2. YOUT
%     3. MF_cmpts
%     4. MF_allcmpts
%     5. difftime
%     6. elapsed_time

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG 
global BDELTA SDELTA SEQ OGSEPER dPGSE
  
bvalues = experiment.bvalues;
qvalues = experiment.qvalues;
gdir = experiment.gdir;
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;
ODEsolve_rtol = experiment.rtol;
ODEsolve_atol = experiment.atol;
yes = 1;  no = 0;

disp(['Magnetization: Solving Bloch-Torrey PDE']);

UG = gdir';
UG = UG/norm(UG);

if (mymesh.Ncmpt == 1 | abs(max(kappa_bdys)) <= 1e-16)
    DO_COUPLING = no;
else
    DO_COUPLING = yes;
end

nexperi = length(sdeltavec);
%disp(['Simulating ',num2str(nexperi),' experiments']);
nb = size(qvalues,2);
elapsed_time = zeros(nb, nexperi);

for icmpt = 1:mymesh.Ncmpt
    %disp(['Working on cmpt ', num2str(icmpt)]);
    coordinates = mymesh.Pts_cmpt_reorder{icmpt};
    elements = mymesh.Ele_cmpt_reorder{icmpt};
    facets = [];
    GX = sqrt(-1)*UG*coordinates;
    FEM_MAT{icmpt}.Q = sparse(length(coordinates),length(coordinates));
    for iboundary = 1:mymesh.Nboundary
        if (kappa_bdys(iboundary) ~= 0)
            neumann = mymesh.Fac_boundary_reorder{icmpt}{iboundary}';
            neumann_nodes = unique(neumann);
            %zvalues = mymesh.Pts_cmpt_reorder{1}(3,neumann_nodes);
            %mvalues = sqrt((20-abs(zvalues)));
            coeffs_flux_matrix = zeros(max(neumann_nodes),1);
            coeffs_flux_matrix(neumann_nodes) = kappa_bdys(iboundary); %*mvalues
            if ~isempty(neumann)
                FEM_MAT{icmpt}.Q = FEM_MAT{icmpt}.Q + flux_matrixP1_3D(neumann,coordinates',coeffs_flux_matrix);
            end
        end
    end
    [FEM_MAT{icmpt}.K,volumes] = stiffness_matrixP1_3D(elements',coordinates',DIFF_cmpts(icmpt));
    FEM_MAT{icmpt}.M = mass_matrixP1_3D(elements',volumes);
    FEM_MAT{icmpt}.A = mass_matrixP1_3D(elements',volumes,GX');
    
    % Add T2 to the BTPDE by Van-Dang Nguyen May 13 2019
    FEM_MAT{icmpt}.K = FEM_MAT{icmpt}.K + 1/T2_cmpts(icmpt)*FEM_MAT{icmpt}.M;

    % calculate the IC
    IC_Pts{icmpt} = IC_cmpts(icmpt)*ones(size(mymesh.Pts_cmpt_reorder{icmpt},2),1);
end

if (DO_COUPLING == yes)
    tic
    [FEMcouple_MAT,FEMcouple_ind0,FEMcouple_indf] ...
        = generate_FEM_coupling(FEM_MAT,mymesh.Ncmpt,mymesh.Nboundary,...
        mymesh.Pts_cmpt_reorder,mymesh.Ele_cmpt_reorder,mymesh.Pts_ind,mymesh.Pts_boundary_reorder,mymesh.Fac_boundary_reorder);
    %stop
    toc
    % IC for coupled case
    IC_couple = zeros(size(FEMcouple_MAT.M,1),1);
    for icmpt = 1:mymesh.Ncmpt
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
    if (SEQ == dPGSE)
        TE = 2*(SDELTA+BDELTA);
    end;
    omega = 2*pi*npervec(iexperi)/SDELTA;
    OGSEPER = 1./omega*2*pi;%% set up number for OGSE
    
    disp(['    Experiment ',num2str(iexperi)]);
    %disp([SDELTA,BDELTA,npervec(iexperi)]);
    
    TLIST = [0,TE];
    
    USE_MIDPOINT = 1;
    
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
                       
            if (USE_MIDPOINT==1)
                options.tol = 1e-6;
                options.dt = 100;
                disp('***Coupled: start mid-point solver'); tic
                sol=midpoint_solver(TLIST,IC_couple,options);
                disp('***Coupled: end mid-point solver'); toc
            else     
                options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                    'Jacobian',@odejac_bt_includeb);            
                disp('***Coupled: start ode solve ode23t'); tic
                sol = ode23t(@odefun_bt_includeb,TLIST,IC_couple,options);
                disp('***Coupled: end ode solve ode23t'); toc
            end;
            
            for icmpt = 1:mymesh.Ncmpt
                YOUT{iexperi}{ib}{icmpt} = sol.y(FEMcouple_ind0(icmpt):FEMcouple_indf(icmpt),:);
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
            end
        else
            %% Solving for case of no coupling between compartments.            
            for icmpt = 1:mymesh.Ncmpt
                FEM_M = FEM_MAT{icmpt}.M;
                FEM_K = FEM_MAT{icmpt}.K;
                FEM_A = FEM_MAT{icmpt}.A;
                FEM_Q = FEM_MAT{icmpt}.Q;
                FEM_G = sparse(zeros(size(FEM_M,1),1));
                ICC = IC_Pts{icmpt};

                if (max(abs(ICC))<=1e-16)
                    sol.y = zeros(size(ICC,1),size(TLIST,2));
                    sol.x = TLIST;
                else
                    if (USE_MIDPOINT==1)
                        options.tol = 1e-6;
                        options.dt = 100;
                        disp('***Uncoupled: start mid-pint solver'); tic
                        sol=midpoint_solver(TLIST,ICC,options);
                        disp('***Uncoupled: end mid-point solver'); toc
                    else
                        options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',ODEsolve_rtol,'Vectorized','on','Stats','off',...
                            'Jacobian',@odejac_bt_includeb);
                        disp('***Uncoupled: start ode solver ode23t'); tic
                        sol = ode23s(@odefun_bt_includeb,TLIST,ICC,options);                
                        disp('***Uncoupled: end ode solver ode23t'); toc
                    end
                end;
                YOUT{iexperi}{ib}{icmpt} = sol.y;
                TOUT{iexperi}{ib}{icmpt} = sol.x;
                MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
%                 if (1==2) 
%                     tp = TLIST(1); % previous time
%                     dt = 100;
%                     t = tp;
%                     u = IC_Pts{icmpt}; 
%                     tol = 1e-6; maxit = 10000;
%                     u_tarray = [u];
%                     tarray = [t];
%                     while t < TLIST(2) + dt
%                         Al = FEM_M - dt/2*(-FEM_K - QVAL*seqprofile(t) *FEM_A);
%                         Ar = FEM_M + dt/2*(-FEM_K - QVAL*seqprofile(tp)*FEM_A);
%                         [u,flag,relres,iter,resvec] = bicgstab(Al,Ar*u, tol, maxit, [], [], u);
%                         u_tarray = [u_tarray, u];
%                         tp = t;
%                         t = t + dt;
%                         tarray = [tarray, t];
%                     end;
%                     YOUT{iexperi}{ib}{icmpt} = u_tarray;
%                     TOUT{iexperi}{ib}{icmpt} = tarray;
%                     MT{iexperi}{ib}{icmpt} = sum(FEM_MAT{icmpt}.M*YOUT{iexperi}{ib}{icmpt},1);
%                     
%                     ICC = IC_Pts{icmpt};
%                     options.tol = 1e-6;
%                     sol=midpoint_solver(TLIST,ICC,options);
%                     e=max(abs(sol.y - u_tarray));
%                     if e>1e-16
%                         stop;
%                     end;
%                 end;
            end;            
        end
        elapsed_time(ib, iexperi)=etime(clock, b_start_time);
    end
end

for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);  
    nb = length(bvec);
    for ib = 1:nb
        for icmpt = 1:mymesh.Ncmpt
            MF_cmpts(icmpt,iexperi,ib) = MT{iexperi}{ib}{icmpt}(end);
            M0(icmpt,iexperi,ib) = MT{iexperi}{ib}{icmpt}(1);
        end
        MF_allcmpts(iexperi,ib) = 0;
        for icmpt = 1:mymesh.Ncmpt
            MF_allcmpts(iexperi,ib) = MF_allcmpts(iexperi,ib) + MF_cmpts(icmpt,iexperi,ib);
        end
        M0_allcmpts(iexperi,ib) = 0;
        for icmpt = 1:mymesh.Ncmpt
            M0_allcmpts(iexperi,ib) = M0_allcmpts(iexperi,ib) + M0(icmpt,iexperi,ib);
        end
    end   
    ib0 = find(abs(bvec)<=1e-16);
    ibn0 = find(abs(bvec)>1e-16);
    if (length(ib0) >= 1)
        for icmpt = 1:mymesh.Ncmpt
            S0(icmpt,iexperi) = MF_cmpts(icmpt,iexperi,ib0(1));
        end
        S0_allcmpts(iexperi) = MF_allcmpts(iexperi,ib0(1));
    else
        S0(1:mymesh.Ncmpt,iexperi) = nan;
        S0_allcmpts(iexperi) = nan;
    end
end