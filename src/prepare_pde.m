function params_domain = prepare_pde(params_cells, params_domain)
%PREPARE_PDE Create initial data for solving PDE.
%
%   params_cells: struct
%  	params_domain: struct
%
%   params_domain: struct


% Extract compartment information
ncell = params_cells.ncell;
shape = params_cells.shape;
include_in = params_cells.include_in;
in_ratio = params_cells.in_ratio;
include_ecs = params_cells.ecs_shape ~= "no_ecs";
ecs_ratio = params_cells.ecs_ratio;

% Check for correct radius ratios and that neurons do not have in-compartments
assert(~include_in || 0 < in_ratio && in_ratio < 1 && shape ~= "neuron");
assert(~include_ecs || 0 < ecs_ratio);

% Number of compartments
ncompartment = (1 + include_in) * ncell + include_ecs;

% Find number of boundaries
switch shape
    case "cylinder"
        % An axon has a side interface, and a top-bottom boundary
        nboundary = (2 * include_in + 1 + include_ecs) * ncell + include_ecs;
    case "sphere"
        % For a sphere, there is one interface
        nboundary = (include_in + 1) * ncell + include_ecs;
    case "neuron"
        % For a neuron, there is one interface
        nboundary = 1 + include_ecs;
end

boundary_markers = false(ncompartment, nboundary);
sz = size(boundary_markers);

if include_in
    % Add in-compartments and in-out-interfaces
    compartments = repmat("in", 1, ncell);
    boundaries = repmat("in,out", 1, ncell);
    boundary_markers(sub2ind(sz, 1:2*ncell, repmat(1:ncell, 1, 2))) = true;
    ncompartment_old = ncell;
    nboundary_old = ncell;
else
    compartments = [];
    boundaries = [];
    ncompartment_old = 0;
    nboundary_old = 0;
end

% Add out-compartments
compartments = [compartments repmat("out", 1, ncell)];

if include_ecs
    % Add ecs-compartment and out-ecs interfaces
    compartments = [compartments "ecs"];
    boundaries = [boundaries repmat("out,ecs", 1, ncell)];
    boundary_markers(sub2ind(sz, ...
        [ncompartment_old+1:ncompartment_old+ncell, repelem(ncompartment, 1, ncell)], ...
        repmat(nboundary_old+1:nboundary_old+ncell, 1, 2))) = true;
    nboundary_old = nboundary_old + ncell;
end

if shape == "cylinder"
    if include_in
        % Add in boundary
        boundaries = [boundaries repmat("in", 1, ncell)];
        ncompartment_old = ncell;
        boundary_markers(sub2ind(sz, 1:ncell, nboundary_old+1:nboundary_old+ncell)) = true;
        nboundary_old = nboundary_old + ncell;
    else
        ncompartment_old = 0;
    end
    
    % Add out boundary
    boundaries = [boundaries repmat("out", 1, ncell)];
    boundary_markers(sub2ind(sz, ncompartment_old+1:ncompartment_old+ncell, ...
        nboundary_old+1:nboundary_old+ncell)) = true;

    if include_ecs
        % Add ecs boundary
        boundaries = [boundaries "ecs"];
        boundary_markers(end, end) = true;
    end
elseif include_ecs
    % Add ecs boundary
    boundaries = [boundaries "ecs"];
    boundary_markers(end, end) = true;
else
    % Add out boundary
    boundaries = [boundaries "out"];
    boundary_markers(end, end) = true;
end

% Diffusion coefficients (tensorize if scalars)
params_domain.diffusivity_in = params_domain.diffusivity_in * eye(3);
params_domain.diffusivity_out = params_domain.diffusivity_out * eye(3);
params_domain.diffusivity_ecs = params_domain.diffusivity_ecs * eye(3);
diffusivity = zeros(3, 3, ncompartment);
diffusivity(:, :, compartments == "in") = repmat(params_domain.diffusivity_in, 1, 1, sum(compartments == "in"));
diffusivity(:, :, compartments == "out") = repmat(params_domain.diffusivity_out, 1, 1, sum(compartments == "out"));
diffusivity(:, :, compartments == "ecs") = repmat(params_domain.diffusivity_ecs, 1, 1, sum(compartments == "ecs"));

% T2-relaxation coefficients
relaxation = zeros(1, ncompartment);
relaxation(compartments == "in") = params_domain.relaxation_in;
relaxation(compartments == "out") = params_domain.relaxation_out;
relaxation(compartments == "ecs") = params_domain.relaxation_ecs;

% Initial conditions
initial_density = zeros(1, ncompartment);
initial_density(compartments == "in") = params_domain.initial_density_in;
initial_density(compartments == "out") = params_domain.initial_density_out;
initial_density(compartments == "ecs") = params_domain.initial_density_ecs;

% Permeability coefficients
permeability = zeros(1, nboundary);
permeability(boundaries == "in") = params_domain.permeability_in;
permeability(boundaries == "out") = params_domain.permeability_out;
permeability(boundaries == "ecs") = params_domain.permeability_ecs;
permeability(boundaries == "in,out") = params_domain.permeability_in_out;
permeability(boundaries == "out,ecs") = params_domain.permeability_out_ecs;


% Update domain parameters with new variables
params_domain.diffusivity = diffusivity;
params_domain.relaxation = relaxation;
params_domain.initial_density = initial_density;
params_domain.permeability = permeability;
params_domain.compartments = compartments;
params_domain.boundaries = boundaries;
params_domain.boundary_markers = boundary_markers;
