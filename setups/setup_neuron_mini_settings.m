%% File name to load or store cell description, surface geometry, mesh, and simulation results
setup.name = "mesh_files/spindle/whole_neurons/16o_spindle13aFI.stl";

%% Geometry parameters
setup.geometry.cell_shape = "neuron";                   % Cell shape; "sphere", "cylinder" or "neuron"
setup.geometry.ncell = 1;                               % Number of cells

%% PDE parameters
setup.pde.diffusivity_out = 0.002;                      % Diffusion coefficient OUT (scalar or 3x3-tensor)

%% Gradient sequences
setup.gradient.values = 10000;                          % g-, q-, or b-values [1 x namplitude]
setup.gradient.values_type = "b";                       % Type of values: "g", "q", or "b" or "g"
setup.gradient.sequences = {                            % Gradient sequences {1 x nsequence}
    PGSE(5000, 10000)
};
setup.gradient.directions = unitsemicircle(1);          % Gradient directions [3 x ndirection]

%% BTPDE experiment parameters (comment block to skip experiment)
setup.btpde = struct;                                   % define an empty setup.btpde to trigger btpde

%% BTPDE midpoint experiment parameters (comment block to skip experiment)
setup.btpde_midpoint.implicitness = 0.5;                % Theta-parameter: 0.5 for Crank-Nicolson
setup.btpde_midpoint.timestep = 2;                      % Time step dt

%% HADC experiment parameters (comment block to skip experiment)
setup.hadc = struct;                                    % define an empty setup.hadc to trigger hadc

%% MF experiment parameters (comment block to skip experiment)
setup.mf.neig_max = Inf;                                % Requested number of eigenvalues
setup.mf.length_scale = 1;                              % Minimum length scale of eigenfunctions

%% Karger model parameters (comment block to skip experiment)
setup.karger.ndirection = 50;                           % Number of directions to compute diffusion tensor
setup.karger.reltol = 1e-4;                             % Relative tolerance for ODE solver
setup.karger.abstol = 1e-6;                             % Absolute tolerance for ODE solver

%% Analytical experiment parameters (comment block to skip experiment)
setup.analytical.length_scale = 0.3;                    % Minimum length scale of eigenfunctions
setup.analytical.eigstep = 1e-8;                        % Minimum distance between eigenvalues
