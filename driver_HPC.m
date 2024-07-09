%DRIVER_HPC A driver that is suitable for massive simulations running in a HPC cluster.
% Add SpinDoctor
clear
restoredefaultpath
addpath(genpath("src"));

% Define fixed inputs
% Geometry parameters
setup.geometry.cell_shape = "neuron";                   % Cell shape; "sphere", "cylinder" or "neuron"
setup.geometry.ncell = 1;                               % Number of cells
setup.geometry.ecs_shape = "tight_wrap";                % Shape of ECS: "no_ecs", "box", "convex_hull", or "tight_wrap".
setup.geometry.ecs_ratio = 0.1;                         % ECS gap (times rmean)
setup.geometry.refinement = -1;                         % Tetgen refinement parameter (comment for automatic)
% PDE parameters
setup.pde.diffusivity_out = 0.002;                      % Diffusion coefficient OUT (scalar or 3x3-tensor)
setup.pde.diffusivity_ecs = 0.002;                      % Diffusion coefficient ECS (scalar or 3x3-tensor)
setup.pde.initial_density_out = 1.0;                    % Initial density in OUT
setup.pde.initial_density_ecs = 1.0;                    % Initial density in ECS
setup.pde.permeability_out = 0;                         % Permeability OUT boundary
setup.pde.permeability_ecs = 0;                         % Permeability ECS boundary
% Gradient sequences
setup.gradient.values = [1000, 3000, 5000, 10000];      % g-, q-, or b-values [1 x namplitude]
setup.gradient.values_type = "b";                       % Type of values: "g", "q", or "b" or "g"
setup.gradient.sequences = {                            % Gradient sequences {1 x nsequence}
    PGSE(12900, 21800)
};
setup.gradient.directions = unitsemisphere(5);         % Gradient directions [3 x ndirection]
% HADC experiment parameters (comment block to skip experiment)
setup.hadc.ode_solver = @ode15s;                        % ODE solver for HADC
setup.hadc.reltol = 1e-4;                               % Relative tolerance for ODE solver
setup.hadc.abstol = 1e-4;                               % Absolute tolerance for ODE solver
% BTPDE experiment parameters (comment block to skip experiment)
setup.btpde.ode_solver = @ode15s;                       % ODE solver for BTPDE
setup.btpde.reltol = 1e-4;                              % Relative tolerance for ODE solver
setup.btpde.abstol = 1e-6;                              % Absolute tolerance for ODE solver
% MF experiment parameters (comment block to skip experiment)
setup.mf.neig_max = Inf;                               % Requested number of eigenvalues
setup.mf.length_scale = 1;                              % Minimum length scale of eigenfunctions
setup.mf.ninterval = 1000;                               % Number of intervals to discretize time profile in MF (if not PGSE and doublePGSE)
setup.mf.eigs.tolerance = 1e-10;                             % Convergence tolerance of eigs
setup.mf.eigs.maxiter = 1000;                                % Maximum number of eigs iterations
% other fixed settings
magnetization_flag = false;
saveroot = 'new_save_path';

% Define variables
name_list = [
    "mesh_files/spindle/whole_neurons/16o_spindle13aFI.stl"
    "mesh_files/spindle/separated_neurons/03a_spindle2aFI_dendrites_1.stl"
    "mesh_files/pyramidal/whole_neurons/03b_pyramidal2aACC.stl"
    "mesh_files/pyramidal/whole_neurons/02a_pyramidal2aFI"
];
kappa_out_ecs = [0, 0.5, 1, 2, 3, 4]*1e-4;
relax_out = [Inf, 10000, 8000, 5000];
relax_ecs = [Inf, 20000, 16000, 10000];
allinds = [length(name_list), length(kappa_out_ecs), length(relax_out)];
results = cell(allinds);

% run simulations
for iall = 1:prod(allinds)
    [iname, ikappa, irelax] = ind2sub(allinds, iall);

    setup.name = name_list(iname);
    setup.pde.permeability_out_ecs = kappa_out_ecs(ikappa);
    setup.pde.relaxation_out = relax_out(irelax);
    setup.pde.relaxation_ecs = relax_ecs(irelax);

    % run_simulation(setup, magnetization_flag, saveroot); % Save results in saveroot
    results{iname, ikappa, irelax} = run_simulation(setup, magnetization_flag, saveroot);
end

%% Postprocess results
% save('results.mat', 'results', '-v7.3');
