function [results, setup, femesh] = load_simulation(setup, magnetization_flag, saveroot)
%LOAD_SIMULATION Load simulation defined by steup

if nargin == 1
    magnetization_flag = true;
    saveroot = 'saved_simul';
elseif nargin == 2
    saveroot = 'saved_simul';
end

% Prepare simulation
[setup, femesh, ~]  = prepare_simulation(setup);

% Solve BTPDE
if isfield(setup, "btpde")
    disp("Loading the BTPDE signals");
    savepath = create_savepath(setup, "btpde", saveroot);
    results.btpde = load_btpde(setup, savepath, magnetization_flag);
end


%% Solve BTPDE using midpoint method
if isfield(setup, "btpde_midpoint")
    disp("Loading the BTPDE midpoint signals");
    savepath = create_savepath(setup, "btpde_midpoint", saveroot);
    results.btpde_midpoint = load_btpde_midpoint( ...
        setup, savepath, magnetization_flag);
end


%% Solve HADC model
if isfield(setup, "hadc")
    disp("Loading the homogenized apparent diffusion coefficient");
    savepath = create_savepath(setup, "hadc", saveroot);
    results.hadc = load_hadc(femesh, setup, savepath);
end


%% Laplace eigendecomposition
if isfield(setup, "mf")
    % Laplace eigendecomposition
    eigenpath = create_savepath(setup, "lap_eig", saveroot);
    results.lap_eig = load_laplace_eig(eigenpath, setup.mf, setup.pde.mean_diffusivity);

    % Compute MF magnetization and signal
    savepath = create_savepath(setup, "mf", saveroot);
    results.mf = load_mf(setup, savepath, magnetization_flag);

    results.mf_hadc = solve_mf_hadc(femesh, setup, results.lap_eig);
end
end
