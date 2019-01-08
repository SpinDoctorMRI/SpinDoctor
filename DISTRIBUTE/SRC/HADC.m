function [ADC_DE,ADC_DE_allcmpts,elapsed_time] ...
    = HADC(experiment,mymesh,DIFF_cmpts,IC_cmpts)
	
	% diffusion equation (zero IC) to get the time-dependent diffusion coefficient


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
                coordVERTICES=zeros(size(neumann,1),3,3);
                coordVERTICES(1:size(neumann,1),1:3,1)=coordinates_t(neumann(:,1),:);
                coordVERTICES(1:size(neumann,1),1:3,2)=coordinates_t(neumann(:,2),:);
                coordVERTICES(1:size(neumann,1),1:3,3)=coordinates_t(neumann(:,3),:);
                [coordNORMALSnew,coordVERTICESnew] = COMPUTE_mesh_normals(coordVERTICES);
				
				ALL_VERTICES=zeros(size(elements,2),3,4);
                
				ALL_VERTICES(1:size(elements,2),1:3,1)=coordinates_t(elements(1,:),:);
                ALL_VERTICES(1:size(elements,2),1:3,2)=coordinates_t(elements(2,:),:);
                ALL_VERTICES(1:size(elements,2),1:3,3)=coordinates_t(elements(3,:),:);        
                ALL_VERTICES(1:size(elements,2),1:3,4)=coordinates_t(elements(4,:),:);        

                midpoint_facet_coords = mean(coordVERTICES(:,:,:),3);
                midpoint_ele_coords = mean(ALL_VERTICES(:,:,:),3);
                for ifacet=1:size(midpoint_facet_coords,1)
                  % find the ending point of the normal of each facet
                  endpoint_of_normal = midpoint_facet_coords(ifacet,:) + coordNORMALSnew(ifacet,:);
                  % compare the position of the ending point of the normal to the
                  % facet
                  endpoint_eval =   coordNORMALSnew(ifacet,1)*(endpoint_of_normal(1)-midpoint_facet_coords(ifacet,1))+...
                                    coordNORMALSnew(ifacet,2)*(endpoint_of_normal(2)-midpoint_facet_coords(ifacet,2))+...
                                    coordNORMALSnew(ifacet,3)*(endpoint_of_normal(3)-midpoint_facet_coords(ifacet,3));

                  % find the closest element to the ending point of the normal  
                  d=sum((midpoint_ele_coords-ones(size(midpoint_ele_coords,1),1)*midpoint_facet_coords(ifacet,:)).^2,2);
                  [mval,mid]=min(d);
                  closestpoint_to_normal = midpoint_ele_coords(mid,:);

                  % compare the position of the closest point to the facet
                  closestpoint_eval = coordNORMALSnew(ifacet,1)*(closestpoint_to_normal(1)-midpoint_facet_coords(ifacet,1))+...
                                      coordNORMALSnew(ifacet,2)*(closestpoint_to_normal(2)-midpoint_facet_coords(ifacet,2))+...
                                      coordNORMALSnew(ifacet,3)*(closestpoint_to_normal(3)-midpoint_facet_coords(ifacet,3));
                  normal_orientation=(endpoint_eval*closestpoint_eval);
                  coordNORMALSnew(ifacet,:) =sign(-normal_orientation)*coordNORMALSnew(ifacet,:);
                end;
                if (1==0)
                    figure;
                    PLOT_3D_stl_patch(coordVERTICESnew,-coordNORMALSnew);
                  
                    MD = mean(coordVERTICESnew,3);
                    sqrt(sum((MD+0.1*coordNORMALSnew).^2,2))-sqrt(sum((MD+0.0*coordNORMALSnew).^2,2))
                    view(3)
                end;
                mycoeff = (coordNORMALSnew(:,1)*UG(1) + coordNORMALSnew(:,2)*UG(2)+coordNORMALSnew(:,3)*UG(3));
                GG = flux_matrixP1_3D(neumann,coordinates', DIFF_cmpts(icmpt)*mycoeff);     
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