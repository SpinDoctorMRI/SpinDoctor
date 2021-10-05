function pde = prepare_pde(setup)
%PREPARE_PDE Create initial data for solving PDE.
%
%   setup: struct
%
%   pde: struct


% Extract compartment information
ncell = setup.geometry.ncell;
cell_shape = setup.geometry.cell_shape;
include_in = setup.geometry.include_in;
in_ratio = setup.geometry.in_ratio;
include_ecs = setup.geometry.ecs_shape ~= "no_ecs";
ecs_ratio = setup.geometry.ecs_ratio;
pde = setup.pde;

% Check for correct radius ratios and that neurons do not have in-compartments
assert(~include_in || 0 < in_ratio && in_ratio < 1 && cell_shape ~= "neuron");
assert(~include_ecs || 0 < ecs_ratio);

% Number of compartments
ncompartment = (1 + include_in) * ncell + include_ecs;

% Find number of boundaries
switch cell_shape
    case "cylinder"
        % An axon has a side interface, and a top-bottom boundary
        nboundary = (include_in + 1) * 2 * ncell + include_ecs;
    case "sphere"
        % For a sphere, there is one interface
        nboundary = (include_in + 1) * ncell + include_ecs;
    case "neuron"
        % For a neuron, there is one interface
        nboundary = 1 + include_ecs;
end

% the list `boundaries` must have the same order as surfaces.facetmarkers
% cylinder: `boundaries` = [('in,out'), ('out,ecs')/('out'), ('in'), 'out', ('ecs')]
% sphere, neuron: `boundaries` = [('in,out'), ('out,ecs')/('out'), ('ecs')]
compartments = [];
boundaries = [];

if include_in
    % Add in-compartments and in-out-interfaces
    compartments = repmat("in", 1, ncell);
    boundaries = repmat("in,out", 1, ncell);
end

% Add out-compartments
compartments = [compartments repmat("out", 1, ncell)];

if include_ecs
    % Add ecs-compartment and out-ecs interfaces
    compartments = [compartments "ecs"];
    boundaries = [boundaries repmat("out,ecs", 1, ncell)];
else
    % Add outer cylinder side wall or sphere/neuron out boundaries
    boundaries = [boundaries repmat("out", 1, ncell)];
end

if cell_shape == "cylinder"
    if include_in
        % Add inner cylinder top and bottom boundary
        boundaries = [boundaries repmat("in", 1, ncell)];
    end
    % Add outer cylinder top and bottom boundary
    boundaries = [boundaries repmat("out", 1, ncell)];
    % Add ecs boundary
    if include_ecs
        boundaries = [boundaries "ecs"];
    end
end

if ismember(cell_shape, ["sphere", "neuron"]) && include_ecs
    % Add ecs boundary
    boundaries = [boundaries "ecs"];
end

% Diffusion coefficients (tensorize if scalars)
pde.diffusivity_in = pde.diffusivity_in * eye(3);
pde.diffusivity_out = pde.diffusivity_out * eye(3);
pde.diffusivity_ecs = pde.diffusivity_ecs * eye(3);
diffusivity = zeros(3, 3, ncompartment);
diffusivity(:, :, compartments == "in") = repmat(pde.diffusivity_in, 1, 1, sum(compartments == "in"));
diffusivity(:, :, compartments == "out") = repmat(pde.diffusivity_out, 1, 1, sum(compartments == "out"));
diffusivity(:, :, compartments == "ecs") = repmat(pde.diffusivity_ecs, 1, 1, sum(compartments == "ecs"));

% T2-relaxation coefficients
relaxation = zeros(1, ncompartment);
relaxation(compartments == "in") = pde.relaxation_in;
relaxation(compartments == "out") = pde.relaxation_out;
relaxation(compartments == "ecs") = pde.relaxation_ecs;

% Initial conditions
initial_density = zeros(1, ncompartment);
initial_density(compartments == "in") = pde.initial_density_in;
initial_density(compartments == "out") = pde.initial_density_out;
initial_density(compartments == "ecs") = pde.initial_density_ecs;

% Permeability coefficients
permeability = zeros(1, nboundary);
permeability(boundaries == "in") = pde.permeability_in;
permeability(boundaries == "out") = pde.permeability_out;
permeability(boundaries == "ecs") = pde.permeability_ecs;
permeability(boundaries == "in,out") = pde.permeability_in_out;
permeability(boundaries == "out,ecs") = pde.permeability_out_ecs;


% Update domain parameters with new variables
pde.diffusivity = diffusivity;
pde.relaxation = relaxation;
pde.initial_density = initial_density;
pde.permeability = permeability;
pde.compartments = compartments;
pde.boundaries = boundaries;
