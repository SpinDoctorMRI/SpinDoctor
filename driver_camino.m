%DRIVER_SPINDOCTOR Solve BTPDE or MF for a sequence inputted from a camino file.
%   Compare different ADC. Plot results in many directions.
%
%   It is highly recommended to read this driver to understand the workflow
%   of SpinDoctor.
%
%   The user is advised to read the latest version
%   from \url{https://github.com/jingrebeccali/SpinDoctor}

% Add SpinDoctor
addpath(genpath("src"));


%% Define inputs

% Get setup
addpath(genpath("setups"));
setup_camino;

if ~isdir(save_path)
    mkdir(save_path)
end

%% Prepare simulation
[setup, femesh, ~, ~]  = prepare_simulation(setup);

%% Perform small experiments
% Short time approximation (STA) of the ADC
% Free diffusion signal
free = compute_free_diffusion(setup.gradient.bvalues, setup.pde.diffusivity, ...
    femesh.volumes, setup.pde.initial_density);

tic
% Perform BTPDE experiments
if isfield(setup, "btpde")
    % Solve BTPDE
    btpde = solve_btpde(femesh, setup,save_path,true);
end
toc
% Perform MF experiments
% Perform Laplace eigendecomposition
tic
if isfield(setup,"mf")
lap_eig = compute_laplace_eig(femesh, setup.pde, setup.mf,save_path);
    
% Compute MF magnetization
mf = solve_mf_cpu(femesh, setup, lap_eig);
end
toc

