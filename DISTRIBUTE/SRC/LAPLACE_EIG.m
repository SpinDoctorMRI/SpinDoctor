function [EIG_value,EIG_proj,EIG_func,ctime] ...
    = LAPLACE_EIG(Pts_cmpt_reorder,Ele_cmpt_reorder,...
    DIFF,VOL,EigLim,EigLim_min,EigIntervals)

% computes the Laplace eigenvalues, eigenfunctions 
% and the first order moments of the products of the eigenfuntctions
% 
% Input:
%     1. Pts_cmpt_reorder
%     2. Ele_cmpt_reorder
%     3. DIFF
%     4. VOL
%     5. EigLim
%     6. EigLim_min
%     7. EigIntervals
%     
% Output:
%     1. EIG_value
%     2. EIG_proj
%     3. EIG_func
%     4. ctime

global FEM_M FEM_K

ttt=cputime;

model_orig = createpde();
geometryFromMesh(model_orig,Pts_cmpt_reorder,Ele_cmpt_reorder);

sigma0 = DIFF;
% compute the eigenvalues and eigenfunctions and save them
specifyCoefficients(model_orig,'m',0,'d',1,'c',sigma0,'a', 0,'f',0);
nface = model_orig.Geometry.NumFaces; %number of face
%     normal dot (c grad u) + q*u = g
for iface = 1:nface
    applyBoundaryCondition(model_orig,'neumann','face',iface,'g',0,'q',0,'Vectorized','on'); % 3-D geometry
end

model_FEM_matrices = assembleFEMatrices(model_orig);
% Define coefficients of PDE (diffusion equation)
FEM_M = model_FEM_matrices.M;
FEM_K = model_FEM_matrices.K;

P = model_orig.Mesh.Nodes;

CenterMass_cmpts = zeros(3,1);
CenterMass_cmpts(1,1) = sum(FEM_M*((P(1,:)')))/VOL;
CenterMass_cmpts(2,1) = sum(FEM_M*((P(2,:)')))/VOL;
CenterMass_cmpts(3,1) = sum(FEM_M*((P(3,:)')))/VOL;

%disp(['center of mass = ']);
%CenterMass_cmpts{icmpt}

P(1,:) = P(1,:)-CenterMass_cmpts(1,1);
P(2,:) = P(2,:)-CenterMass_cmpts(2,1);
P(3,:) = P(3,:)-CenterMass_cmpts(3,1);


model_shift = createpde();
geometryFromMesh(model_shift,P,Ele_cmpt_reorder);
sigma0 = DIFF;
% compute the eigenvalues and eigenfunctions and save them
specifyCoefficients(model_shift,'m',0,'d',1,'c',sigma0,'a', 0,'f',0);
nface = model_shift.Geometry.NumFaces; %number of face
%     normal dot (c grad u) + q*u = g
for iface = 1:nface
    applyBoundaryCondition(model_shift,'neumann','face',iface,'g',0,'q',0,'Vectorized','on'); % 3-D geometry
end
D = 1;
C = sigma0;
A = 0;

nint = EigIntervals;


V = [];
eigvals = [];
evec = linspace(0,EigLim,nint+1);
evec(1) = EigLim_min;
for iint = 2:nint+1
    e0 = evec(iint-1)+1e-16;
    ef = evec(iint);
    disp([e0,ef]/DIFF);
    EigenResults = solvepdeeig(model_shift,[e0,ef])
    V = [V,EigenResults.Eigenvectors];
    eigvals = [eigvals; EigenResults.Eigenvalues];
end


neig = length(eigvals);

eigvals = eigvals(1:neig);

V = V(:,1:neig);

neig = length(eigvals);

% normalize the eigenfunctions and calculate
Ueig_N = V;
for ieig = 1:neig
    normUeig(ieig) = sqrt(Ueig_N(:,ieig)'*FEM_M*Ueig_N(:,ieig));
    Ueig_N(:,ieig) = Ueig_N(:,ieig)/normUeig(ieig);
end

EIG_proj = zeros(3,neig,neig);
for ieig = 1:neig
    for jeig = 1:neig
        EIG_proj(1,ieig,jeig) = sum(FEM_M*(Ueig_N(:,jeig).*Ueig_N(:,ieig).*(P(1,:)')));
        EIG_proj(2,ieig,jeig) = sum(FEM_M*(Ueig_N(:,jeig).*Ueig_N(:,ieig).*(P(2,:)')));
        EIG_proj(3,ieig,jeig) = sum(FEM_M*(Ueig_N(:,jeig).*Ueig_N(:,ieig).*(P(3,:)')));
    end
end


% normalized eigenfunctions
EIG_func = Ueig_N;
% eigenvalues
EIG_value = eigvals;

ctime = cputime-ttt;


