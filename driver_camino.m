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
[setup, femesh, surfaces, cells]  = prepare_simulation(setup);

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

%% Plots
do_plots = true;

if ~do_plots
    return
end

% Plot surface triangulation
plot_surface_triangulation(surfaces);

% Plot the finite element mesh
plot_femesh(femesh, setup.pde.compartments);
% plot_femesh_everywhere(femesh, "");

% Plot information about the geometry
plot_geometry_info(setup, femesh);


% Plot BTPDE magnetization in some directions
if isfield(setup, "btpde")

    for iseq = 1:setup.nsequence
         field = btpde.magnetization(iseq);
         title_str = sprintf(...
                    "BTPDE magnetization. Sequence %d of %d", ...
                    iseq, setup.nsequence);
         plot_field_everywhere(femesh, field, title_str);
    end
    clear field
end

if isfield(setup, "mf")
    % Relative error between BTPDE and MF signal
    signal_allcmpts_relerr = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ...
        ./ max(abs(btpde.signal_allcmpts), [], 3);

    % Difference between BTPDE and MF signal, normalized by initial signal
    signal_allcmpts_abserr_vol = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ./ setup.pde.initial_signal;
    
    % Display difference
    disp("Error (normalized by initial signal)");
    disp(signal_allcmpts_abserr_vol);
    disp("Relative error:");
    disp(signal_allcmpts_relerr);

end
