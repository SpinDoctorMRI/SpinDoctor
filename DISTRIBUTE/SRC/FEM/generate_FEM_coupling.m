function [FEMcouple_MAT,ind0,indf] ...
    = generate_FEM_coupling_JRL(FEM_MAT,Ncmpt,Nboundary,...
    Pts_cmpt_reorder,Ele_cmpt_reorder,Pts_ind,Pts_boundary_reorder,Fac_boundary_reorder)
  
  for iboundary = 1:Nboundary
    for icmpt = 1:Ncmpt
      ord{icmpt}{iboundary} = (Pts_ind{icmpt}(Pts_boundary_reorder{icmpt}{iboundary}));
      Npts_bdy(icmpt,iboundary) = length(ord{icmpt}{iboundary});
    end
  end
  
  ndof_total = 0;
  ndof = zeros(Ncmpt,1);
  ind0 = zeros(Ncmpt,1);
  indf = zeros(Ncmpt,1);
  for icmpt = 1:Ncmpt
    ndof(icmpt) = size(FEM_MAT{icmpt}.M,1);%number of unknown in one cylinder
    ndof_total = ndof_total+ndof(icmpt); % total number of unknown in domain
  end
  FEMcouple_MAT.M = sparse((ndof_total),(ndof_total));
  FEMcouple_MAT.K = sparse((ndof_total),(ndof_total));
  FEMcouple_MAT.A = sparse((ndof_total),(ndof_total));
  FEMcouple_MAT.Q = sparse((ndof_total),(ndof_total));
  for icmpt = 1:Ncmpt % take the FE matrix into right place 
    ind0(icmpt) = sum(ndof(1:icmpt-1))+1;
    indf(icmpt) = sum(ndof(1:icmpt));
    FEMcouple_MAT.M(ind0(icmpt):indf(icmpt),ind0(icmpt):indf(icmpt)) = FEM_MAT{icmpt}.M;
    FEMcouple_MAT.K(ind0(icmpt):indf(icmpt),ind0(icmpt):indf(icmpt)) = FEM_MAT{icmpt}.K;      
    FEMcouple_MAT.A(ind0(icmpt):indf(icmpt),ind0(icmpt):indf(icmpt)) = FEM_MAT{icmpt}.A;
    FEMcouple_MAT.Q(ind0(icmpt):indf(icmpt),ind0(icmpt):indf(icmpt)) = FEM_MAT{icmpt}.Q;
  end  
  for iboundary = 1:Nboundary
    cmpts_touch = find(Npts_bdy(:,iboundary) ~= 0);
    ntmp = length(cmpts_touch);
    if (ntmp > 2) 
      disp('each interface touch only 1 or 2 cmpts');
      stop
    elseif (ntmp == 2)
      cmpt1 = cmpts_touch(1);
      cmpt2 = cmpts_touch(2);
          
      Q1 = FEM_MAT{cmpt1}.Q;
      Q2 = FEM_MAT{cmpt2}.Q;
      Q12 = zeros(size(Q1,1),size(Q2,2));
      Q21 = zeros(size(Q2,1),size(Q1,2));
      N = Npts_bdy(cmpt1,iboundary);
      for ii = 1:N
        pt_ind_orig = ord{cmpt1}{iboundary}(ii);
        jj = find(ord{cmpt2}{iboundary}==pt_ind_orig);
        ind2 = Pts_boundary_reorder{cmpt2}{iboundary}(jj);
        ind1 = Pts_boundary_reorder{cmpt1}{iboundary}(ii);
        Q12(:,ind2) = -Q1(:,ind1);
        Q21(:,ind1) = -Q2(:,ind2);
      end  
      FEMcouple_MAT.Q(ind0(cmpt1):indf(cmpt1),ind0(cmpt2):indf(cmpt2)) = sparse(Q12);
      FEMcouple_MAT.Q(ind0(cmpt2):indf(cmpt2),ind0(cmpt1):indf(cmpt1)) = sparse(Q21);      
      %disp(['interface boundary: ',num2str(iboundary),', connecting cmpts ',num2str([cmpt1,cmpt2])]);
    elseif (ntmp == 1)
      % zero out the boundary condition on boundary that is not interface between two compartments.
      cmpt1 = cmpts_touch(1);
      N = Npts_bdy(cmpt1,iboundary);
      icmpt = cmpt1;    
      for ii = 1:N
        ind1 = Pts_boundary_reorder{cmpt1}{iboundary}(ii);      
        FEMcouple_MAT.Q(ind0(icmpt):indf(icmpt),ind0(icmpt)-1+ind1) = 0;
      end    
      %disp(['out boundary: ',num2str(iboundary),', only cmpt ',num2str([cmpt1])]);
    end
  end
    
