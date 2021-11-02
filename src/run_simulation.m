function results = run_simulation(setup, magnetization_flag, saveroot)
%RUN_SIMULATION Run simulation defined by steup

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
    disp("Computing or loading the BTPDE signals");
    savepath = create_savepath(setup, "btpde", saveroot);
    results.btpde = solve_btpde(femesh, setup, savepath, magnetization_flag);
end


%% Solve BTPDE using midpoint method
if isfield(setup, "btpde_midpoint")
    disp("Computing or loading the BTPDE midpoint signals");
    savepath = create_savepath(setup, "btpde_midpoint", saveroot);
    results.btpde_midpoint = solve_btpde_midpoint( ...
        femesh, setup, savepath, magnetization_flag ...
    );
end


%% Solve HADC model
if isfield(setup, "hadc")
    disp("Computing or loading the homogenized apparent diffusion coefficient");
    savepath = create_savepath(setup, "hadc", saveroot);
    results.hadc = solve_hadc(femesh, setup, savepath);
end


%% Laplace eigendecomposition
if isfield(setup, "mf")
    % Laplace eigendecomposition
    eigenpath = create_savepath(setup, "lap_eig", saveroot);
    results.lap_eig = compute_laplace_eig(femesh, setup.pde, setup.mf, eigenpath);

    % Compute MF magnetization and signal
    savepath = create_savepath(setup, "mf", saveroot);
    results.mf = solve_mf(femesh, setup, lap_eig, savepath, magnetization_flag);

    results.mf_hadc = solve_mf_hadc(femesh, setup, lap_eig);
end

% Solve Karger model
if isfield(setup, "karger")
    % Solve analytical analytical model
    results.karger = solve_karger(femesh, setup);
end

% Perform analytical experiment
if isfield(setup, "analytical")
    % Solve analytical analytical model
%     volumes = compute_cmptwise_volume(femesh.volumes, setup.pde.compartments);
%     analytical_signal = solve_analytical(setup, volumes); % With FE volumes
    results.analytical_signal = solve_analytical(setup); % With pure volumes
end
end
