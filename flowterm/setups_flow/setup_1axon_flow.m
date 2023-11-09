
%% File name to load or store cell description, surface geometry, mesh, and simulation results
setup.name = "flowterm/mesh_files/cylinders/1axon";

%% Geometry parameters
setup.geometry.cell_shape = "cylinder";                 % Cell shape: "sphere", "cylinder" or "neuron"
setup.geometry.ncell = 1;                               % Number of cells
setup.geometry.rmin = 2.1;                                % Minimum radius of cells
setup.geometry.rmax = 2.1;                                % Maximum radius of cells
setup.geometry.dmin = 0.1;                              % Minimum distance between cells (times mean(rmin,rmax))
setup.geometry.dmax = 0.2;                              % Maximum distance between cells (times mean(rmin,rmax))
setup.geometry.height = 100;%2*75;%20;                             % Cylinder height (ignored if not cylinder)
setup.geometry.deformation = [0.0; 0.0];            % Domain deformation; [a_bend, a_twist]

%% Geometry parameters
setup.geometry.include_in = false;                      % Ratio Rin/R, within range [0,0.99]
setup.geometry.in_ratio = 0.1;                          % Ratio Rin/R, within range [0,0.99]
setup.geometry.ecs_shape = "no_ecs";                % Shape of ECS: "no_ecs", "box", "convex_hull", or "tight_wrap".
setup.geometry.ecs_ratio = 0.1;                         % ECS gap; percentage in side length

%% Finite element mesh parameters (comment block to use default parameters)
setup.geometry.refinement  = 1;%-1;%0.1;                     % Tetgen refinement parameter (comment for automatic)

%% PDE parameters
setup.pde.diffusivity_in = 0e-3;                       % Diffusion coefficient IN (scalar or 3x3-tensor)
setup.pde.diffusivity_out = 2e-3;                      % Diffusion coefficient OUT (scalar or 3x3-tensor)
setup.pde.diffusivity_ecs = 0e-3;                      % Diffusion coefficient ECS (scalar or 3x3-tensor)
setup.pde.relaxation_in = Inf;                          % T2-relaxation IN. No relaxation: Inf
setup.pde.relaxation_out = Inf;                         % T2-relaxation OUT. No relaxation: Inf
setup.pde.relaxation_ecs = Inf;                         % T2-relaxtion ECS. No relaxation: Inf
setup.pde.initial_density_in = 1.0;                     % Initial density in IN
setup.pde.initial_density_out = 1.0;                    % Initial density in OUT
setup.pde.initial_density_ecs = 1.0;                    % Initial density in ECS
setup.pde.permeability_in_out = 0;                   % Permeability IN-OUT interface
setup.pde.permeability_out_ecs = 0;                  % Permeability OUT-ECS interface
setup.pde.permeability_in = 0;                          % Permeability IN boundary
setup.pde.permeability_out = 0;                         % Permeability OUT boundary
setup.pde.permeability_ecs = 0;                         % Permeability ECS boundary

%% Gradient sequences
setup.gradient.values = [200];%[0];                    % g-, q-, or b-values [1 x namplitude]
setup.gradient.values_type = "g";                       % Type of values: "g", "q", or "b" or "g"
setup.gradient.sequences = {                            % Gradient sequences {1 x nsequence} 
PGSE(1e3, 1e3)
}';
setup.gradient.directions = [0.0; 0.0; 1.0];% Gradient directions [3 x ndirection]

%% BTPDE midpoint experiment parameters (comment block to skip experiment)
setup.btpde_midpoint.implicitness = 0.5;              % Theta-parameter: 0.5 for Crank-Nicolson
setup.btpde_midpoint.timestep = 0.5;                    % Time step dt