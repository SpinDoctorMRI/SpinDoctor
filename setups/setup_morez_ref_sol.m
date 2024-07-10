%SETUP_NEURON Define setup structure for SpinDoctor.
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
%% Geometry parameters
setup.geometry.cell_shape = "neuron";                   % Cell shape; "sphere", "cylinder" or "neuron"
setup.geometry.ncell = 1;                               % Number of cells
setup.geometry.deformation = [0; 0];                    % Domain deformation; [a_bend,a_twist]
setup.geometry.include_in = false;                      % Ratio Rin/R, within range [0,0.99]
setup.geometry.in_ratio = 0.7;                          % Ratio Rin/R, within range [0,0.99]
setup.geometry.ecs_shape = "no_ecs";                    % Shape of ECS: "no_ecs", "box", "convex_hull", or "tight_wrap".
setup.geometry.ecs_ratio = 0.2;                         % ECS gap; percentage in side length

% setup.geometry.refinement = 10;                       % Tetgen refinement parameter (comment for automatic)
setup.geometry.tetgen_options = "-pq1.2aVCn";              % Tetgen options (priority is inferior to refinement)

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
setup.pde.permeability_out_ecs = 1e-4;                  % Permeability OUT-ECS interface
setup.pde.permeability_in = 0;                          % Permeability IN boundary
setup.pde.permeability_out = 0;                         % Permeability OUT boundary
setup.pde.permeability_ecs = 0;                         % Permeability ECS boundary

%% Gradient sequences
seq = read_scheme('custom_sequence_files/morez_Q3_camino.scheme');
% if isfile('custom_sequence_files/morez_Q3.txt')
% b_tensors = read_b_tensors('custom_sequence_files/morez_Q3.txt');
% for iseq = 1:length(seq)
%     seq{iseq}.b_tensor = b_tensors(:,iseq);
% end
% end
setup.gradient.sequences =seq(2:11);
%% BTPDE experiment parameters (comment block to skip experiment)
setup.btpde.ode_solver = @ode15s;                       % ODE solver for BTPDE
setup.btpde.reltol = 1e-4;                              % Relative tolerance for ODE solver
setup.btpde.abstol = 1e-6;                              % Absolute tolerance for ODE solver
setup.btpde.rerun = false;    
%% MF experiment parameters (comment block to skip experiment)
% Length scale hard-coded for these experiments from the diffusivity values and sequence length.
char_length_scale = sqrt(2*3*setup.pde.diffusivity_in*1000*101);
setup.mf.length_scale =char_length_scale/5;             % Minimum length scale of eigenfunctions
setup.mf.neig_max = 2000;                               % Requested number of eigenvalues
setup.mf.ninterval = 200;                               % Number of intervals to discretize time profile in MF (if not PGSE and doublePGSE)
setup.mf.eigs.sigma = 1e-8;    
setup.mf.rerun=false;                                        % Rerun simulation with or without saved results
