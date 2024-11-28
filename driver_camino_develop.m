%DRIVER_CAMINO Solve BTPDE or MF for a sequence inputted from a camino file.
%
%   Provides a template to demonstrate how camino scheme sequences can be
%   integrated into SpinDoctor
%
%   The user is advised to read the latest version
%   from \url{https://github.com/jingrebeccali/SpinDoctor}

% Add SpinDoctor
addpath(genpath("src"));
addpath(genpath("drivers_postprocess"));

%% Define inputs
mesh =  "C:\Users\amcsween\SpinDoctor_mesh_files/spindle/whole_neurons/04b_spindle3aFI.ply";

%% Run simulation
[results,femesh,~,~]= run_simulations_neuron(mesh,"setup_camino_develop");

mf = results.mf_cell;
setup = results.setup; nsequence = setup.nsequence;
disp("Camino sequence results are stored in mf");

%% Plot Geometry
do_plots = true;

if do_plots
    
    % Plot the finite element mesh
    plot_femesh(femesh, setup.pde.compartments);
    % plot_femesh_everywhere(femesh, "");
    
    % Plot information about the geometry
    plot_geometry_info(setup, femesh);
end


%% Plot Sequences
% Plot sequence
for iseq = 1%:nsequence
    title_str = sprintf("Sequence %d of %d",iseq,setup.nsequence);
    setup.gradient.sequences{iseq}.plot(title_str);
end

gamma = setup.gamma * 1e-6;

% Study B_tensor
for iseq=1%:nsequence
    B_tensor = setup.gradient.sequences{iseq}.get_B_tensor(gamma);
    fprintf("For sequence %d of %d:\n",iseq,setup.nsequence);
    B_tensor
    plot_tensor(B_tensor,sprintf("B tensor for sequence %d of %d",iseq, setup.nsequence));
end


%% Fitting Diffusion tensor
% See Q-space trajectory imaging for multidimensional diffusion MRI of the
% human brain, (2016) Westin

% Note that for a fixed sequence, B is a rank 2 tensor and we also compute
% $B^(\otimes 2)$, which is rank 4. However, we represent it by a rank two
% tensor as is done in Westin 2016. This is done by represening B by a
% vector (so rank 1) and then $B^(\otimes 2)$ by the outer product of this
% vector. The representations are chosen to preserve the inner product of
% tensors.

%% Extract B-tensor data and signals

% Map from symmetric rank 2 tensor to vector representation (see eq 7)
B_vector = @(B) [B(1,1) B(2,2) B(3,3) sqrt(2)*B(1,2) sqrt(2)*B(1,3) sqrt(2)*B(2,3)];

% Matrix of weights to represent rank 4 tensor by rank 2 tensor with
% reduced degrees of freedom (see eq 11)
W= [
    1 1 1 sqrt(2) sqrt(2) sqrt(2);
    1 1 1 sqrt(2) sqrt(2) sqrt(2);
    1 1 1 sqrt(2) sqrt(2) sqrt(2);
    sqrt(2) sqrt(2) sqrt(2) 2 2 2;
    sqrt(2) sqrt(2) sqrt(2) 2 2 2;
    sqrt(2) sqrt(2) sqrt(2) 2 2 2
    ] ;
% Indices of unique entries in $B^(\otimes 2)$.
ind = [1 2 3 4 5 6 8 9 10 11 12  15 16 17 18 22 23 24 29 30 36];


B_tensors = zeros(3,3,nsequence);
b_vectors = zeros(6,nsequence);
B2_tensors = zeros(21,nsequence);
for iseq=1:nsequence
    B= setup.gradient.sequences{iseq}.get_B_tensor(gamma);
    B_tensors(:,:,iseq) = B;
    b = B_vector(B);
    B2 = b'*b;
    symm = W.*B2;
    B2_tensors(:,iseq) = reshape(symm(ind),[],1);
    b_vectors(:,iseq) = b;
end
B_tensors = reshape(B_tensors,9,120);

S = log(real(mf.signal)./femesh.total_volume)';



%% (Diffusion tensor) (EQ 12) Solve with Matlab
% Set input variables
X = b_vectors';

% Solve equation, representing D,D^(\otimes 2) as a vector of size 30
A = X\S;

% Split A into D and D^(\otimes 2) vectors
d = -A(1:6);

% Diffusion matrix from vector representation
Diffusion_matrix = @(d) [d(1) d(6)/sqrt(2) d(5)/sqrt(2); d(6)/sqrt(2) d(2) d(4)/sqrt(2); d(5)/sqrt(2) d(4)/sqrt(2) d(3)];
D = Diffusion_matrix(d);
rel_signal_error = abs(exp(X*A) - exp(S))./abs(exp(S));
abs_signal_error = abs(exp(X*A) - exp(S));


disp("Diffusion tensor");
disp("Diffusion matrix = ");
D
disp("Eigenvalues = ");
eig(D)

fprintf("Mean relative signal error = %f, Std dev relative signal error = %f\n" + ...
    "Mean absolute signal error = %f, Std dev absolute signal error = %f\n", ...
    mean(rel_signal_error),std(rel_signal_error),mean(abs_signal_error),std(abs_signal_error));


%% (Diffusion tensor and covariance) (EQ 13) Solve with Matlab
% Set input variables
X = [b_vectors',B2_tensors'];

% Solve equation, representing D,D^(\otimes 2) as a vector of size 30
A = X\S;

% Split A into D and D^(\otimes 2) vectors
d = -A(1:6);
c = 2*A(7:27);

% Diffusion matrix from vector representation
Diffusion_matrix = @(d) [d(1) d(6)/sqrt(2) d(5)/sqrt(2); d(6)/sqrt(2) d(2) d(4)/sqrt(2); d(5)/sqrt(2) d(4)/sqrt(2) d(3)];
D = Diffusion_matrix(d);
C = Covariance_matrix(c);

rel_signal_error = abs(exp(X*A) - exp(S))./abs(exp(S));
abs_signal_error = abs(exp(X*A) - exp(S));

disp("Diffusion tensor with covariance tensor")
disp("Diffusion tensor = ")
D
disp("Eigenvalues = ");
eig(D)
disp("Diffusion covariance tensor = ")
C
disp("Eigenvalues = ");
eig(C)

fprintf("Mean relative signal error = %f, Std dev relative signal error = %f\n" + ...
    "Mean absolute signal error = %f, Std dev absolute signal error = %f\n", ...
    mean(rel_signal_error),std(rel_signal_error),mean(abs_signal_error),std(abs_signal_error));

%% Misc variables for below
G1 = [1 0 0 0 0 0; 0 0 0 0 0 1/sqrt(2); 0 0 0 0 1/sqrt(2) 0];
G2 = [0 0 0 0 0 1/sqrt(2); 0 1 0 0 0 0; 0 0 0 1/sqrt(2) 0 0];
G3 = [0 0 0 0 1/sqrt(2) 0; 0 0 0 1/sqrt(2) 0 0; 0 0 1 0 0 0];

disp("User needs to set path to cvx here.")
cd ../cvx
cvx_setup
cd ../SpinDoctor

%% (Diffusion tensor) (EQ 12) Solve with cvx

X = b_vectors';
S = log(real(mf.signal)./femesh.total_volume)';
clear d_cvx;
cvx_begin
    variable d_cvx(6);
    minimize(norm( S - X*d_cvx));
    subject to
        [G1*d_cvx G2*d_cvx G3*d_cvx] == semidefinite(3);
cvx_end

D_cvx = Diffusion_matrix(d_cvx);
disp("Diffusion tensor with cvx");
disp("Diffusion matrix = ");
D_cvx
disp("Eigenvalues = ");
eig(D_cvx)

rel_signal_error = abs(exp(X*d_cvx) - exp(S))./abs(exp(S));
abs_signal_error = abs(exp(X*d_cvx) - exp(S));

fprintf("Mean relative signal error = %f, Std dev relative signal error = %f\n" + ...
    "Mean absolute signal error = %f, Std dev absolute signal error = %f\n", ...
    mean(rel_signal_error),std(rel_signal_error),mean(abs_signal_error),std(abs_signal_error));

%% (Diffusion tensor and covariance) (EQ 13) Solve with cvx

X = [b_vectors',B2_tensors'];
S = log(real(mf.signal)./femesh.total_volume)';
clear d_cvx c_cvx;
cvx_begin
    variable d_cvx(6);
    variable c_cvx(21);
    minimize(norm( S - X*[d_cvx' c_cvx']'));
    subject to
        [G1*d_cvx G2*d_cvx G3*d_cvx] == semidefinite(3);
        Covariance_matrix(c_cvx) == semidefinite(6);
cvx_end

D_cvx = Diffusion_matrix(d_cvx);
C_cvx = Covariance_matrix(c_cvx)
disp("Diffusion tensor with cvx");
disp("Diffusion matrix = ");
D_cvx
disp("Eigenvalues = ");
eig(D_cvx)
disp("Diffusion covariance tensor = ")
C_cvx
disp("Eigenvalues = ");
eig(C_cvx)
rel_signal_error = abs(exp(X*[d_cvx' c_cvx']') - exp(S))./abs(exp(S));
abs_signal_error = abs(exp(X*[d_cvx' c_cvx']') - exp(S));


fprintf("Mean relative signal error = %f, Std dev relative signal error = %f\n" + ...
    "Mean absolute signal error = %f, Std dev absolute signal error = %f\n", ...
    mean(rel_signal_error),std(rel_signal_error),mean(abs_signal_error),std(abs_signal_error));

%% Compute biological parameters
% Note I use Dmean2 to indicate $< D >^{\otimes 2}$ and 
%            D2mean to indicate $< D^{\otimes 2} >$

% Need to ensure that bounds are met by model.

C = C_cvx;

% Equations 16-19
E_iso  = 1/3 * eye(6);
E_bulk = 1/9*[ones(3,3) zeros(3,3);zeros(3,3) zeros(3,3)];
E_shear = E_iso - E_bulk;
V_md = dot(C(:),E_bulk(:));

% Equation 6
Dmean2 =squeeze(tensorprod(d,d));

D2mean = C + Dmean2;

% Equation 28 
C_mu = 1.5* dot(D2mean(:),E_shear(:))/dot(D2mean(:),E_iso(:));
% Equation 29 microscopic anistropy
mu_FA = sqrt(C_mu);

% Equation 30
C_M = 1.5* dot(Dmean2(:),E_shear(:))/dot(Dmean2(:),E_iso(:));
% Equation 30 macroscopic anistropy
FA = sqrt(C_M);

% Equation 32 orientation parameter
OP = sqrt( dot(Dmean2(:),E_shear(:))/dot(D2mean(:),E_shear(:)));

% Equation 33 microscopic orientation coherence
C_c = C_M/C_mu;

% Equation 41 bulk kurtosis
K_bulk = 3*dot(C(:),E_bulk(:))/dot(Dmean2(:),E_bulk(:));

% Equation 42 shear kurtosis
K_shear= (6/5)*dot(C(:),E_shear(:))/dot(Dmean2(:),E_bulk(:));

% Equation 43 microscopic kurtosis
K_mu= (6/5)*dot(D2mean(:),E_shear(:))/dot(Dmean2(:),E_bulk(:));


fprintf("microscopic anisotropy = %f\n", mu_FA);
fprintf("macroscopic anisotropy = %f\n", FA);
fprintf("orientation parameter = %f\n",OP );
fprintf("bulk kurtosis = %f\n", K_bulk);
fprintf("shear kurtosis = %f\n", K_shear);
fprintf("microscopic kurtosis = %f\n", K_mu);

fprintf("Bounds model should have:\n")
fprintf("microscopic anisotropy in [0,1]\n");
fprintf("orientation parameter in [0,1]\n");
fprintf("May also be more, double check paper\n")

%% Functions

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
































