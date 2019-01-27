function [ADC_DE,ADC_DE_allcmpts,elapsed_time] ...
    = HADC(experiment,mymesh,DIFF_cmpts,IC_cmpts)

% diffusion equation (zero IC) to get the time-dependent diffusion coefficient
% 
% Input:
%     1. experiment is a structure with 8 elements:
%         ngdir_total 
%         gdir        
%         sdeltavec   
%         bdeltavec   
%         seqvec      
%         npervec    
%         rtol       
%         atol            
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
%     4. IC_cmpts
% 
% Output:
%     1. ADC_DE
%     2. ADC_DE_allcmpts
%     3. elapsed_time

gdir = experiment.gdir;
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;
rtol = experiment.rtol;
atol = experiment.atol;

global FEM_M FEM_K FEM_A FEM_Q FEM_G
global QVAL UG
global BDELTA SDELTA SEQ OGSEPER

disp(['H-ADC Model: solving DE']);

UG = gdir';
UG = UG/norm(UG);
        
for icmpt = 1:mymesh.Ncmpt
    [VOL(icmpt)] ...
        = get_volume_mesh(mymesh.Pts_cmpt_reorder{icmpt},mymesh.Ele_cmpt_reorder{icmpt});
end

VOL_allcmpts = 0;

for icmpt = 1:mymesh.Ncmpt
    VOL_allcmpts  = VOL_allcmpts + VOL(icmpt);
end

for icmpt = 1:mymesh.Ncmpt
    VOL_frac(icmpt) = VOL(icmpt)/VOL_allcmpts;
end

nexperi = length(sdeltavec);
elapsed_time=zeros(mymesh.Ncmpt, nexperi);


for icmpt = 1:mymesh.Ncmpt
    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % To replace the pde matrices
        coordinates = mymesh.Pts_cmpt_reorder{icmpt}; 
        elements = mymesh.Ele_cmpt_reorder{icmpt};    
        [FEM_K,volumes]=stiffness_matrixP1_3D(elements',coordinates',DIFF_cmpts(icmpt));
        FEM_M=mass_matrixP1_3D(elements',volumes);   
        FEM_G=sparse(size(FEM_K,1),1);
        FEM_A=sparse(size(FEM_K,1),size(FEM_K,2));
        FEM_Q=sparse(size(FEM_K,1),size(FEM_K,2)); 
        coordinates_t=coordinates';
        one = sparse(size(FEM_K,1),1);
        one(:,:) = 1;

        for iboundary = 1:mymesh.Nboundary
            neumann = mymesh.Fac_boundary_reorder{icmpt}{iboundary}';
            if sum(size(neumann))>0
                [FacA,FacC,FacN] = get_surfacenormal_mesh(mymesh.Pts_cmpt_reorder{icmpt},mymesh.Ele_cmpt_reorder{icmpt},neumann');
                mycoeff = (FacN(1,:)*UG(1) + FacN(2,:)*UG(2)+FacN(3,:)*UG(3));
                GG = flux_matrixP1_3D(neumann,coordinates', DIFF_cmpts(icmpt)*mycoeff');     
                FEM_G = FEM_G + GG*one;
            end;
        end;
        
        % MyM = model_FEM_matrices{icmpt}.M;
        % MyK = model_FEM_matrices{icmpt}.K;
        % MyA = model_FEM_matrices{icmpt}.A;
        % MyQ = model_FEM_matrices{icmpt}.Q;
        % MyG = DIFF_cmpts(icmpt)*model_FEM_matrices{icmpt}.G;
        % errorG=sum((MyG - FEM_G).^2);
        % errorK=sum(sum((MyK - FEM_K).^2));
        % errorM=sum(sum((MyM - FEM_M).^2));
        % if errorG+errorK+errorM>1e-7
            % disp('Matrices are different!'); 
            % disp(['Error M:', num2str(errorM),', Error K:', num2str(errorK), ', Error G:', num2str(errorG)]);
            % % [full(MyG),  full(FEM_G)]
            % stop;
        % end;
        % disp(['Error M:', num2str(errorM),', Error K:', num2str(errorK), ', Error G:', num2str(errorG)]);
        % To replace the pde matrices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    
    ODEsolve_atol = atol;
    ODEsolve_rtol = rtol;
    
    options = odeset('Mass',FEM_M,'AbsTol',ODEsolve_atol,'RelTol',...
        ODEsolve_rtol,'Vectorized','on','Stats','off',...
        'Jacobian',@odejac_bt_includeb);
    disp(['    Compartment ',num2str(icmpt)]);
    nexperi = length(sdeltavec);  % Is it correct to be here???
    for iexperi = 1:nexperi
		disp(['      Experiment ',num2str(iexperi)]);
        iex_start_time = clock;
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
        
        ADC_DE(icmpt,iexperi) = DIFF_cmpts(icmpt)-a;
        elapsed_time(icmpt, iexperi) = etime(clock, iex_start_time);
    end
end
ADC_DE_allcmpts = nan*ones(nexperi,1);
for iexperi = 1:nexperi
    ADC_DE_allcmpts(iexperi,1) = sum((IC_cmpts.*VOL_frac)'.*ADC_DE(:,iexperi))./sum((IC_cmpts.*VOL_frac)');
end