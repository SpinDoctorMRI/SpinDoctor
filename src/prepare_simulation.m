function [setup, femesh, surfaces, cells] = prepare_simulation(setup)
%PREPARE_SETUP Prepare setup.
%   The parameters are added to or removed from the input structure.
%
%   setup: struct
%
%   setup: struct
%   femesh: struct
%   surfaces: struct
%   cells: struct


% We here assume we are working with water protons
% https://en.wikipedia.org/wiki/Proton_magnetic_moment#:~:text=The%20proton's%20gyromagnetic%20ratio
setup.gamma = 0.2675222005; % unit: rad / (us*mT)

% Check and clean up the geometry entry
setup.geometry = prepare_geometry(setup);

% Set up the PDE model in the geometrical compartments.
setup.pde = prepare_pde(setup);

% Check consistency of gradient sequences and simulation methods
setup = prepare_experiments(setup);

% Update setup info
setup.ncompartment = length(setup.pde.compartments);
setup.nboundary = length(setup.pde.boundaries);

setup.nsequence = length(setup.gradient.sequences);
if isfield(setup.gradient,'directions') && isfield(setup.gradient,'values')
    setup.ndirection = size(setup.gradient.directions, 2);
    setup.namplitude = length(setup.gradient.values);
else
    setup.ndirection = 0;
    setup.namplitude = 0;
end
% Create or load finite element mesh
[femesh, surfaces, cells] = create_geometry(setup);

% Compute initial total signal
setup.pde.initial_signal = setup.pde.initial_density * femesh.volumes';
% Compute mean diffusivity
setup.pde.mean_diffusivity = compute_mean_diffusivity(setup.pde.diffusivity, femesh);
