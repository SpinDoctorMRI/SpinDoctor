function signal = solve_analytical(setup, volumes)
%SOLVE_ANALYTICAL Compute signal for analytical cylinders or spheres.
%   The signal is computed using an analytical matrix formalism exploiting
%   geometrical invariances symmetries, reducing the Laplace eigenvalue problem
%   into finding radial eigenvalues for different angular eigenmodes.
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
% This function converts SpinDoctor units to SI units.
% Units in analytical module - all SI-units:
%   permeability: m / s
%   water diffusion coefficient: m^2 / s
%   sphere radius: m
%   gyromagnetic ratio: rad / T / s
%   gradient pulse duration: s
%   time between pulses: s
%
%   setup: struct
%   volumes: [1 x ncompartment]
%
%   signal: [namplitude x nsequence x ndirection]


% Extract parameters
qvalues = setup.gradient.qvalues;
sequences = setup.gradient.sequences;
rmean = (setup.geometry.rmin + setup.geometry.rmax) / 2;

% Sizes
namplitude = length(setup.gradient.values);
nsequence = length(setup.gradient.sequences);

if setup.geometry.cell_shape == "cylinder"
    dim = 2;
    get_vol = @(r) pi * r.^2 * setup.geometry.height;
else
    dim = 3;
    get_vol = @(r) 4 * pi / 3 * r.^3;
end

% Get OUT parameters
rho_out = setup.pde.initial_density_out;
r_out = (setup.geometry.rmin + setup.geometry.rmax) / 2;
D_out = trace(setup.pde.diffusivity_out) / 3;
T_out = setup.pde.relaxation_out;
W_out = 0;

% Get IN parameters
if setup.geometry.include_in
    rho_in = setup.pde.initial_density_in;
    r_in = setup.geometry.in_ratio * r_out;
    D_in = trace(setup.pde.diffusivity_in) / 3;
    W_in = setup.pde.permeability_in_out;
    T_in = setup.pde.relaxation_in;
else
    % Empty IN
    rho_in = [];
    r_in = [];
    D_in = [];
    W_in = [];
    T_in = [];
end

% Get ECS parameters
if setup.geometry.ecs_shape ~= "no_ecs"
    % Include ECS
    rho_ecs = setup.pde.initial_density_ecs;
    r_ecs = r_out + setup.geometry.ecs_ratio * rmean;
    D_ecs = trace(setup.pde.diffusivity_ecs) / 3;
    W_ecs = setup.pde.permeability_ecs;
    T_ecs = setup.pde.relaxation_ecs;

    % Add interface permeability between out and ecs
    W_out = setup.pde.permeability_out_ecs;
else
    % Empty ECS
    rho_ecs = [];
    r_ecs = [];
    D_ecs = [];
    W_ecs = [];
    T_ecs = [];
    
    if dim == 3
        % Add surface relaxivity for the outermost sphere ("out")
        W_out = setup.pde.permeability_out;
    end
end

% Create parameters structure
d = dim;
r = [r_in r_out r_ecs] * 1e-6;
D = [D_in D_out D_ecs] * 1e-6;
W = [W_in W_out W_ecs] * 1;
T = [T_in T_out T_ecs] * 1e-6;
rho = [rho_in rho_out rho_ecs];
L = length(D);

params.d = d;
params.r = r;
params.D = D;
params.W = W;
params.T = T;

% Get exact volumes
volumes_exact = get_vol(params.r * 1e6);

% Subtract inner volumes (each volume is contained within the next)
volumes_exact(2:end) = volumes_exact(2:end) - volumes_exact(1:end-1);


% Volumes
if nargin < nargin(@solve_analytical)
    disp("Using exact volumes")
    volumes = volumes_exact;
else
    disp("Exact volumes:")
    disp(volumes_exact)
    disp("Finite element volumes:")
    disp(volumes)
    disp("Relative error:")
    disp(abs(volumes - volumes_exact) ./ volumes_exact)
    disp("Using finite element volumes")
end

% Volume weighted quantities
D_mean = volumes / sum(volumes) * params.D';
S0 = volumes * rho.';

% Upper Laplace eigenvalue limit
eiglim = pi^2 * D_mean / setup.analytical.length_scale^2;

% Compute radial and angular eigenvalues (lambda and nu), with
% alpha^2 = lambda and n^2 = nu for cylinders and n(n+1) = nu for spheres
[alpha, n] = compute_eigenvalues(sqrt(eiglim) * 1e6, params, sqrt(setup.analytical.eigstep) * 1e3);

N = length(alpha);
lam = alpha.^2;

beta = zeros(1, N);
for m = 1:N
    I = compute_I(alpha(m), alpha(m), n(m), params);
    beta(m) = 1 / sqrt(2 * sum(I));
end

J = zeros(N, 1);
for m = 1:N
    if n(m) == 0
        J(m) = compute_J(alpha(m), params);
    end
end

% Eigenfunction integrals
U = (2 * (d == 2) + sqrt(6) * (d == 3)) * J .* beta';

% Diffusion matrix
Lam = diag(lam);

% Moment matrix
B = zeros(N, N);
for m1 = 1:N
    for m2 = m1:N
        if abs(n(m1) - n(m2)) == 1
            K = compute_K(alpha(m1), n(m1), alpha(m2), n(m2), params);
            if d == 2
                eps = sqrt(1 + (n(m1) == 0) + (n(m2) == 0));
            else
                eps = (n(m1) + n(m2) + 1) / sqrt((2 * n(m1) + 1) * (2 * n(m2) + 1));
            end
            B(m1, m2) = eps * beta(m1) * beta(m2) * K;
            B(m2, m1) = B(m1, m2);
        end
    end
end

% Mass matrices
Bri = zeros(N, N, L);
for m1 = 1:N
    for m2 = m1:N
        if n(m1) == n(m2)
            I = compute_I(alpha(m1), alpha(m2), n(m1), params);
            Bri(m1, m2, :) = 2 * beta(m1) * beta(m2) * I;
            Bri(m2, m1, :) = Bri(m1, m2, :);
        end
    end
end

% Relaxation matrix
Br = zeros(N, N);
for i = 1:L
    Br = Br + Bri(:, :, i) / params.T(i);
end

% Coefficients of initial density
U(:) = 0;
U(1) = 1;

% Iterate over experiments
signal = zeros(namplitude, nsequence);
for iseq = 1:nsequence
    for iamp = 1:namplitude

        seq = sequences{iseq};
        w = qvalues(iamp, iseq) * 1e12 * params.r(end);

        del = seq.delta * 1e-6;
        Del = seq.Delta * 1e-6;

        % Compute signal attenuation
        E = U' ...
            * expm(-del * (Lam + Br + 1i * w * B)) ...
            * expm(-(Del - del) * (Lam + Br)) ...
            * expm(-del * (Lam + Br - 1i * w * B)) ...
            * U;
        signal(iamp, iseq) = E * S0;
    end
end

disp(5);
