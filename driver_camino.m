%DRIVER_CAMINO Solve BTPDE or MF for a sequence inputted from a camino file.
%
%   Provides a template to demonstrate how camino scheme sequences can be
%   integrated into SpinDoctor
%
%   The user is advised to read the latest version
%   from \url{https://github.com/jingrebeccali/SpinDoctor}

% Add SpinDoctor
addpath(genpath("src"));


%% Define inputs

% Get setup
addpath(genpath("setups"));
setup_camino;

[~,cellname,~] = fileparts(setup.name);
save_path = sprintf("saved_simul/%s",cellname);
fprintf('Saving to %s\n',save_path);

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
mf = solve_mf(femesh, setup, lap_eig,save_path,false);
end
toc

disp("Camino sequence results are stored in btpde.camino and mf.camino");
disp("PGSE sequence results are stored in btpde.const and mf.const (const for constant direction vector).");

