function pde = prepare_pde(setup)
%PREPARE_PDE Create initial data for solving PDE.
%
%   setup: struct
%
%   pde: struct


% Extract compartment information
cell_shape = setup.geometry.cell_shape;
ncell = setup.geometry.ncell;
include_in = setup.geometry.include_in;
include_ecs = setup.geometry.ecs_shape ~= "no_ecs";

pde = setup.pde;

% check diffusivity (tensorize if scalars), relaxation, permeability
pde.diffusivity_out = pde.diffusivity_out * eye(3);
assert(check_diffusion_tensor(pde.diffusivity_out));
if ~isfield(pde, 'initial_density_out')
    pde.initial_density_out = 1.0;
end
if isfield(pde, 'relaxation_out')
    assert(pde.relaxation_out > 0);
else
    pde.relaxation_out = Inf;
end
if isfield(pde, 'permeability_out')
    assert(pde.permeability_out >= 0);
else
    pde.permeability_out = 0;
end

if include_ecs
    pde.diffusivity_ecs = pde.diffusivity_ecs * eye(3);
    assert(check_diffusion_tensor(pde.diffusivity_ecs));
    if ~isfield(pde, 'initial_density_ecs')
        pde.initial_density_ecs = 1.0;
    end
    if isfield(pde, 'relaxation_ecs')
        assert(pde.relaxation_ecs > 0);
    else
        pde.relaxation_ecs = Inf;
    end
    if isfield(pde, 'permeability_ecs')
        assert(pde.permeability_ecs >= 0);
    else
        pde.permeability_ecs = 0;
    end
    if isfield(pde, 'permeability_out_ecs')
        assert(pde.permeability_out_ecs >= 0);
    else
        pde.permeability_out_ecs = 0;
    end
else
    useless_fields = {'diffusivity_ecs', 'initial_density_ecs', ...
        'relaxation_ecs', 'permeability_ecs', 'permeability_out_ecs'};
    pde = rmfields(pde, useless_fields);
end

if include_in
    pde.diffusivity_in = pde.diffusivity_in * eye(3);
    assert(check_diffusion_tensor(pde.diffusivity_in));
    if ~isfield(pde, 'initial_density_in')
        pde.initial_density_in = 1.0;
    end
    if isfield(pde, 'relaxation_in')
        assert(pde.relaxation_in > 0);
    else
        pde.relaxation_in = Inf;
    end
    if isfield(pde, 'permeability_in')
        assert(pde.permeability_in >= 0);
    else
        pde.permeability_in = 0;
    end
    if isfield(pde, 'permeability_in_out')
        assert(pde.permeability_in_out >= 0);
    else
        pde.permeability_in_out = 0;
    end
else
    useless_fields = {'diffusivity_in', 'initial_density_in', ...
        'relaxation_in', 'permeability_in', 'permeability_in_out'};
    pde = rmfields(pde, useless_fields);
end

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

% Initialization
diffusivity = zeros(3, 3, ncompartment);
relaxation = zeros(1, ncompartment);
initial_density = zeros(1, ncompartment);
permeability = zeros(1, nboundary);

% out compartments
diffusivity(:, :, compartments == "out") = repmat(pde.diffusivity_out, 1, 1, sum(compartments == "out"));
initial_density(compartments == "out") = pde.initial_density_out;
relaxation(compartments == "out") = pde.relaxation_out;
permeability(boundaries == "out") = pde.permeability_out;

if include_ecs
    diffusivity(:, :, compartments == "ecs") = repmat(pde.diffusivity_ecs, 1, 1, sum(compartments == "ecs"));
    initial_density(compartments == "ecs") = pde.initial_density_ecs;
    relaxation(compartments == "ecs") = pde.relaxation_ecs;
    permeability(boundaries == "ecs") = pde.permeability_ecs;
    permeability(boundaries == "out,ecs") = pde.permeability_out_ecs;
end

if include_in
    diffusivity(:, :, compartments == "in") = repmat(pde.diffusivity_in, 1, 1, sum(compartments == "in"));
    initial_density(compartments == "in") = pde.initial_density_in;
    relaxation(compartments == "in") = pde.relaxation_in;
    permeability(boundaries == "in") = pde.permeability_in;
    permeability(boundaries == "in,out") = pde.permeability_in_out;
end

% Update domain parameters with new variables
pde.diffusivity = diffusivity;
pde.relaxation = relaxation;
pde.initial_density = initial_density;
pde.permeability = permeability;
pde.compartments = compartments;
pde.boundaries = boundaries;
end

function flag = check_diffusion_tensor(diffusion_tensor)
    % need to check the properties of intrinsic diffusion tensor
    eigvals = eig(diffusion_tensor);
    flag = all(eigvals >= 0) && all(size(diffusion_tensor) == 3);
end
