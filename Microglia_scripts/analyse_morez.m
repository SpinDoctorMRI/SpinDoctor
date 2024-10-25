addpath(genpath('../src'));
mesh= "mesh_files/spindle/04b_spindle3aFI.ply";
tetgen_options="-pq1.2a0.5O9VCn";
addpath(genpath('src'))
addpath(genpath('setups'))
setup_file='setup_morez_microglia';

fid = fopen("cells_human.txt",'r');meshes = textscan(fid,"%s");fclose(fid);
meshes = meshes{1}; meshes = string(meshes);
file_parts = split(meshes,"/");
types = file_parts(:,3);
ncells = length(meshes);
%%

for i = 1%:ncells
    mesh = meshes(i);
    [results,femesh,~,~] = run_simulations_microglia(mesh,setup_file,tetgen_options);
end
mf = results.mf_cell;
sequences = results.setup.gradient.sequences;clc

%%
% Compute B tensors.

W= [
    1 1 1 sqrt(2) sqrt(2) sqrt(2);
    1 1 1 sqrt(2) sqrt(2) sqrt(2);
    1 1 1 sqrt(2) sqrt(2) sqrt(2);
    sqrt(2) sqrt(2) sqrt(2) 2 2 2;
    sqrt(2) sqrt(2) sqrt(2) 2 2 2;
    sqrt(2) sqrt(2) sqrt(2) 2 2 2
    ] ;
% W = W.*(2*triu(ones(6,6)) - eye(6));
nsequences = length(sequences);

B_vector = @(B) [B(1,1) B(2,2) B(3,3) sqrt(2)*B(1,2) sqrt(2)*B(1,3) sqrt(2)*B(2,3)];
bs = zeros(nsequences,6);
b2s = zeros(nsequences,21);
ind = [1 2 3 4 5 6 8 9 10 11 12  15 16 17 18 22 23 24 29 30 36];
gamma = results.setup.gamma * 1e-6;
Bs = zeros(3,3,nsequences);
B2s = zeros(6,6,nsequences);

for iseq =  1:nsequences
seq  = sequences{iseq};
% B_int = seq.get_B_tensor(gamma);
v=  seq.b_tensor(3:end);
axial = squeeze(tensorprod(v',v));
A = [v eye(3)];
A = orth(A);
axis_sym_1 = A(:,2);axis_sym_2 = A(:,3);
axis_sym_1 = squeeze(tensorprod(axis_sym_1',axis_sym_1));
axis_sym_2 = squeeze(tensorprod(axis_sym_2',axis_sym_2));

B= seq.b_tensor(1)*axial +seq.b_tensor(2)*axis_sym_1 +seq.b_tensor(2)*axis_sym_2 ;
b = B_vector(B);
B2 = b'*b;
bs(iseq,:) =b;
% symm = W.*(B2+B2')/2;
symm = W.*B2;
b2s(iseq,:) =reshape(symm(ind),[],1);
% errors(iseq) = norm(B_int - B)/norm(B);
Bs(:,:,iseq) = B;
B2s(:,:,iseq) = B2;
end

B = reshape(Bs,9,120);

B2 = reshape(B2s,36,120);

%% Straight linear regression
% Diffusion_matrix = @(d) [d(1) d(6)/sqrt(2) d(5)/sqrt(2); d(6)/sqrt(2) d(2) d(4)/sqrt(2); d(5)/sqrt(2) d(4)/sqrt(2) d(3)];

% log_S = @(d,c,B,B2) -dot(B(:),reshape(Diffusion_matrix(d),[],1)) + 0.5*dot(B2(:),reshape(Covariance_matrix(c),[],1));
%  X = [bs b2s];
%  % Need to solve AX = B where B are real data and then d = -A(1:6), c = 2*A(7:42)
% Y = log(real(mf.signal)./femesh.total_volume)';
% X = X(2:end,:); Y = Y(2:end);
% 
% A = X \ Y;
% d = -A(1:6);
% c = 2*A(7:42);
% 
% D = Diffusion_matrix(d)
% C = Covariance_matrix(c)
% 
% signal_error = abs(exp(X*A) - exp(Y))./abs(exp(Y));
% 
% positive_definite(A)
%% Straight diffusion tensor imaging ( THis is probably best
% %  X = [bs];
% %  Y = log(real(mf.signal)./femesh.total_volume)';
% % % X = X(2:end,:);Y = Y(2:end);
% % d0 = -X \ Y;
% % D0 = Diffusion_matrix(d0)
% % eig(D0)
% % 
% % C0 = squeeze(tensorprod(d0,d0));
% % 
% % C0 = eye(6);
% % C0 = C0./W;
% % c0 = C0(ind);
clear D_cvx;
cvx_begin
    variable D_cvx(3,3);
    minimize(norm( Y + B'*D_cvx(:)));
    subject to
        D_cvx== semidefinite(3);
cvx_end


error_cvx = norm(Y-B'*D_cvx(:))./norm(Y)

%%
clear D_cvx C_cvx
cvx_begin
    variable D_cvx(3,3);
    variable C_cvx(6,6);
    minimize(norm( Y - (-B'*D_cvx(:) + 1/2 * B2'*C_cvx(:))));
    subject to
        D_cvx == semidefinite(3);
        C_cvx == semidefinite(6);
cvx_end

error_cvx = norm(Y-B'*D_cvx(:))./norm(Y)




%% Straight diffusion tensor imaging (cvx)
X = [bs];
Y = log(real(mf.signal)./femesh.total_volume)';
clear d;
cvx_begin
    variable d(6);
    minimize(norm( Y - X*d));
    subject to
        [G1*d G2*d G3*d] == semidefinite(3);
cvx_end


D = [G1*d G2*d G3*d] 
error_cvx = norm(Y-X*d)./norm(Y)
eig(D)
%%
d0 = -X\Y;
D0 = Diffusion_matrix(d0)
errorerrer= norm(Y-X*d0)./norm(Y)
eig(D0)
%% Matlab optimization
% % % A0  = A;
% % % postive_definite = @(A) [-eig(Diffusion_matrix(A(1:6))) , [] ]
% %  X = [bs b2s];
% %  % Need to solve AX = B where B are real data and then d = -A(1:6), c = 2*A(7:42)
% % Y = log(real(mf.signal)./femesh.total_volume)';
% % % X = X(2:end,:); Y = Y(2:end);
% % % x = fmincon(fun,A0,zeros(size(A))',0,zeros(size(A))',0,0,0,nonlcon)
% % A0 = [d0'  0*c0 ]';
% % positive_definite(A0)
% % 
% % options = optimoptions('fmincon','Display','iter','Algorithm','sqp','StepTolerance',1e-16,'ConstraintTolerance',1e-25, ...
% %     'MaxFunctionEvaluations',1000000,'MaxIterations',40000);
% % problem.options = options;
% % problem.solver = 'fmincon';
% % problem.objective = @(A) norm(X*A- Y);
% % problem.x0 = A0;
% % problem.nonlcon = @positive_definite;
% % [A,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD]= fmincon(problem);
% % d = A(1:6);
% % D = Diffusion_matrix(d)
% % c = A(7:end); C = Covariance_matrix(c);
% % positive_definite(A)
% % OUTPUT.message
% % error = norm(X*A- Y)./norm(Y)

%%
X = [bs b2s];
clear d_cvx c_cvx;
cvx_begin
    variable d_cvx(6);
    variable c_cvx(21);
    minimize(norm( Y - X*[d_cvx' c_cvx']'));
    subject to
        Diffusion_matrix(d_cvx) == semidefinite(3);
        Covariance_matrix(c_cvx) == semidefinite(6);
cvx_end

D_cvx = Diffusion_matrix(d_cvx);
C_cvx = Covariance_matrix(c_cvx);

A_cvx =[d_cvx' c_cvx']';

error_cvx = norm(X*A_cvx- Y)./norm(Y)


%%


S = (-D_cvx(:)'*B + 1/2 * C_cvx(:)'*B2)';
error_cvx = norm(Y -S)/norm(Y)

%%
clear D_cvx C_cvx
cvx_begin
    variable D_cvx(3,3);
    variable C_cvx(6,6);
    minimize(norm( Y - (-D_cvx(:)'*B + 1/2 * C_cvx(:)'*B2)'));
    subject to
        D_cvx == semidefinite(3);
        C_cvx == semidefinite(6);
cvx_end




%%
% Compute microscopic anisotropy.

C = C_cvx;
E_iso  = 1/3 * eye(6);
E_bulk = 1/9*[ones(3,3) zeros(3,3);zeros(3,3) zeros(3,3)];
E_shear = E_iso - E_bulk;

V_md = dot(C(:),E_bulk(:));

D2 = C + squeeze(tensorprod(d,d));

C_mu = 1.5* dot(D2(:),E_shear(:))/dot(D2(:),E_iso(:));

mu_FA = sqrt(C_mu)

D2 =  squeeze(tensorprod(d,d));

C_M = 1.5* dot(D2(:),E_shear(:))/dot(D2(:),E_iso(:));

FA = sqrt(C_M)

C_c = C_M/C_mu

[c,ceq] = positive_definite(A)

%%
% Compute Orientation Parameter.

G1 = [1 0 0 0 0 0; 0 0 0 0 0 1/sqrt(2); 0 0 0 0 1/sqrt(2) 0];
G2 = [0 0 0 0 0 1/sqrt(2); 0 1 0 0 0 0; 0 0 0 1/sqrt(2) 0 0];
G3 = [0 0 0 0 1/sqrt(2) 0; 0 0 0 1/sqrt(2) 0 0; 0 0 1 0 0 0];

% function D = Diffusion_matrix(d)
% D = [d(1) d(6)/sqrt(2) d(5)/sqrt(2); d(6)/sqrt(2) d(2) d(4)/sqrt(2); d(5)/sqrt(2) d(4)/sqrt(2) d(3)];
% 
% end
function C = Covariance_matrix(c)
C = [
    c(1) c(2) c(3) sqrt(2)*c(4) sqrt(2)*c(5) sqrt(2)*c(6);
    c(2) c(7) c(8) sqrt(2)*c(9) sqrt(2)*c(10) sqrt(2)*c(11);
    c(3) c(8) c(12) sqrt(2)*c(13) sqrt(2)*c(14) sqrt(2)*c(15);
    sqrt(2)*c(4) sqrt(2)*c(9) sqrt(2)*c(13) 2*c(16) 2*c(17) 2*c(18);
    sqrt(2)*c(5) sqrt(2)*c(10) sqrt(2)*c(14) 2*c(17) 2*c(19) 2*c(20);
    sqrt(2)*c(6) sqrt(2)*c(11) sqrt(2)*c(15) 2*c(18) 2*c(20) 2*c(21)
];
end


function [c,ceq] = positive_definite(A)
    % c = [-eig(Diffusion_matrix(A(1:6)));-eig(Covariance_matrix(A(7:end)))];
    e1 = [1 0 0]; e2 = [0 1 0]; e3 = [0 0 1];
    f1 = [1 0 0 0 0 0]; f2 = [0 1 0 0 0 0]; f3 = [0 0 1 0 0 0];
    f4 = [0 0 0 1 0 0]; f5 = [0 0 0 0 1 0]; f6 = [0 0 0 0 0 1];
    % N_C_constraints = 100;
    % B = rand(6,N_C_constraints);
    D = Diffusion_matrix(A(1:6));
    C= Covariance_matrix(A(7:end));
    % U = A(7:48);
    % % C = cov_from_vec(U);
    % V = U(6+(1:36));
    
    c = -[e1*D*e1' e2*D*e2' e3*D*e3' 10*f1*C*f1' 10*f2*C*f2' 10*f3*C*f3' 10*f4*C*f4' 10*f5*C*f5' 10*f6*C*f6'];
    % ceq = vecnorm(V,2,1) - 1;
    ceq = [];
end
function e = row(i)
    e= zeros(1,21);
    e(i) = 1;
end