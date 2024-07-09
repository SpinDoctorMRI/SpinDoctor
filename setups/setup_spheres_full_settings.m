%SETUP_SPHERES_FULL_SETTINGS Define setup structure for SpinDoctor.
%
%   The setup structure may contain the following substructures:
%
%       geometry:
%           Geometry parameters. Can contain parameters for specific cell types
%           (spheres, cylinders). Determines whether to include inner and ECS
%           compartments), and finite element mesh parameters.
%   
%       pde:
%           Domain parameters. Contains PDE parameters (material properties)
%           defined for the domains "IN", "OUT", "ECS" and their boundaries
%   
%       gradient:
%           Gradient sequence parameters. It determines the three properties
%           `directions`, `amplitudes` and `sequences`. Available sequences are:
%               PGSE(delta, Delta)
%               DoublePGSE(delta, Delta)
%               CosOGSE(delta, Delta, nperiod)
%               SinOGSE(delta, Delta, nperiod)
%               CustomSequence(delta, Delta, @timeprofile)
%           The directions can be provided manually, or created by the functions
%           `unitsphere`, `unitcircle`, or `unitsemicercle`.
%
%   The precense of any of the following substructures triggers the corresponding experiment:
%
%       btpde:          Solve Bloch-Torrey PDE with P1-FEM and provided ODE
%                       solver
%
%       btpde_midpoint: Solve Bloch-Torrey PDE with P1-FEM and Crank-Nicolson
%                       solver
%
%       hadc:           Solve the equation for the homogenized apparent
%                       diffusion coefficient using P1-FEM
%
%       mf:             Compute the matrix formalism signal
%
%       analytical:     Compute analytical signal one analyticaled sphere or
%                       cylinder using truncated radial matrix formalism
%
%       karger:         Solve for the finite pulse Karger model

%% File name to load or store cell description, surface geometry, mesh, and simulation results
setup.name = "mesh_files/spheres/spheres_full_settings";

%% Geometry parameters
setup.geometry.cell_shape = "sphere";                   % Cell shape; "sphere", "cylinder" or "neuron"
setup.geometry.ncell = 3;                               % Number of cells
setup.geometry.rmin = 1;                                % Minimum radius
setup.geometry.rmax = 3;                                % Maximum radius
setup.geometry.dmin = 0.2;                              % Minimum distance between cells (times mean(rmin,rmax))
setup.geometry.dmax = 0.3;                              % Maximum distance between cells (times mean(rmin,rmax))
setup.geometry.include_in = true;                       % Switch of in-compartment
setup.geometry.in_ratio = 0.3;                          % Ratio Rin/R, within range [0,0.99]
setup.geometry.ecs_shape = "tight_wrap";                % Shape of ECS: "no_ecs", "box", "convex_hull", or "tight_wrap".
setup.geometry.ecs_ratio = 0.2;                         % ECS gap (times rmean)
setup.geometry.refinement = 1;                          % Tetgen refinement parameter (comment for automatic)
% setup.geometry.tetgen_options = "-pAVR";              % Tetgen options (priority is superior to refinement)

%% PDE parameters
setup.pde.diffusivity_in = 0.002;                       % Diffusion coefficient IN (scalar or 3x3-tensor)
setup.pde.diffusivity_out = 0.002;                      % Diffusion coefficient OUT (scalar or 3x3-tensor)
setup.pde.diffusivity_ecs = 0.002;                      % Diffusion coefficient ECS (scalar or 3x3-tensor)
setup.pde.relaxation_in = Inf;                          % T2-relaxation IN. No relaxation: Inf
setup.pde.relaxation_out = Inf;                         % T2-relaxation OUT. No relaxation: Inf
setup.pde.relaxation_ecs = Inf;                         % T2-relaxtion ECS. No relaxation: Inf
setup.pde.initial_density_in = 1.0;                     % Initial density in IN
setup.pde.initial_density_out = 1.0;                    % Initial density in OUT
setup.pde.initial_density_ecs = 1.0;                    % Initial density in ECS
setup.pde.permeability_in_out = 1e-4;                   % Permeability IN-OUT interface
setup.pde.permeability_out_ecs = 3e-4;                  % Permeability OUT-ECS interface
setup.pde.permeability_in = 1e-5;                       % Permeability IN boundary
setup.pde.permeability_out = 1e-5;                      % Permeability OUT boundary
setup.pde.permeability_ecs = 1e-5;                      % Permeability ECS boundary

%% Gradient sequences
setup.gradient.values = 0:2500:10000;                   % g-, q-, or b-values [1 x namplitude]
setup.gradient.values_type = "b";                       % Type of values: "g", "q", or "b" or "g"
setup.gradient.sequences = {                            % Gradient sequences {1 x nsequence}
    PGSE(5000, 10000)
    PGSE(10000, 100000)
};
setup.gradient.directions = unitsemicircle(5);          % Gradient directions [3 x ndirection]

%% BTPDE experiment parameters (comment block to skip experiment)
setup.btpde.ode_solver = @ode15s;                       % ODE solver for BTPDE
setup.btpde.reltol = 1e-4;                              % Relative tolerance for ODE solver
setup.btpde.abstol = 1e-6;                              % Absolute tolerance for ODE solver
setup.btpde.rerun = true;                               % Rerun simulation with or without saved results

%% BTPDE midpoint experiment parameters (comment block to skip experiment)
setup.btpde_midpoint.implicitness = 0.5;                % Theta-parameter: 0.5 for Crank-Nicolson
setup.btpde_midpoint.timestep = 2;                      % Time step dt
setup.btpde.rerun = true;                               % Rerun simulation with or without saved results

%% HADC experiment parameters (comment block to skip experiment)
setup.hadc.ode_solver = @ode15s;                        % ODE solver for HADC
setup.hadc.reltol = 1e-4;                               % Relative tolerance for ODE solver
setup.hadc.abstol = 1e-4;                               % Absolute tolerance for ODE solver
setup.hadc.rerun = true;                                % Rerun simulation with or without saved results

%% MF experiment parameters (comment block to skip experiment)
setup.mf.neig_max = 500;                                % Requested number of eigenvalues
setup.mf.length_scale = 1;                              % Minimum length scale of eigenfunctions
setup.mf.ninterval = 500;                               % Number of intervals to discretize time profile in MF (if not PGSE)
setup.mf.rerun_eigen = true;                            % Rerun eigendecomposition
setup.mf.rerun = true;                                  % Rerun MF simulation
setup.mf.eigs.sigma = -1e-8;                            % Type of eigenvalues (default: all eigs close to -1e-8)
setup.mf.eigs.tolerance = 1e-10;                        % Convergence tolerance of eigs
setup.mf.eigs.maxiter = 1000;                           % Maximum number of eigs iterations
setup.mf.eigs.ssdim = 14000;                            % Maximum size of Krylov subspace

%% Karger model parameters (comment block to skip experiment)
setup.karger.ndirection = 50;                           % Number of directions to compute diffusion tensor
setup.karger.ode_solver = @ode45;                       % ODE solver for BTPDE
setup.karger.reltol = 1e-4;                             % Relative tolerance for ODE solver
setup.karger.abstol = 1e-6;                             % Absolute tolerance for ODE solver

%% Analytical experiment parameters (comment block to skip experiment)
setup.analytical.length_scale = 0.3;                    % Minimum length scale of eigenfunctions
setup.analytical.eigstep = 1e-8;                        % Minimum distance between eigenvalues
