function geometry = prepare_geometry(setup)
%PREPARE_GEOMETRY Prepare or clean up setup.geometry.
%
%   setup: struct
%
%   geometry: struct


% Extract compartment information
geometry = setup.geometry;
cell_shape = geometry.cell_shape;

% check ncell
assert(isnumeric(geometry.ncell) && geometry.ncell >= 1);
% check if geometry.ncell is integer
if geometry.ncell ~= round(geometry.ncell)
    warning("Assign setup.geometry.ncell to %d.", round(geometry.ncell));
    geometry.ncell = round(geometry.ncell);
end

% check cell_shape and geometry parameters
switch cell_shape
    case "cylinder"
        assert(geometry.rmin <= geometry.rmax);
        assert(geometry.dmin <= geometry.dmax);
        assert(0 < geometry.height);

        % check ecs
        geometry = check_ecs(geometry);
        % check in
        geometry = check_in(geometry);
        % check deformation
        geometry = check_deformation(geometry);
    case "sphere"
        assert(geometry.rmin <= geometry.rmax);
        assert(geometry.dmin <= geometry.dmax);
        geometry = rmfields(geometry, {'height'});

        % check ecs
        geometry = check_ecs(geometry);
        % check in
        geometry = check_in(geometry);
        % check deformation
        geometry = check_deformation(geometry);
    case "neuron"
        assert(geometry.ncell == 1);
        useless_fields = {'rmin', 'rmax', 'dmin', 'dmax', ...
            'height', 'in_ratio'};
        geometry = rmfields(geometry, useless_fields);

        % check ecs
        geometry = check_ecs(geometry);
        % check in
        geometry.include_in = false;
        % check deformation
        if isfield(geometry, 'deformation') && any(geometry.deformation)
            warning("Deformation is not available for neurons. Set deformation to [0; 0].")
        end
        geometry.deformation=[0, 0];
    otherwise
        error('Do not support cell shape: %s.', cell_shape);
end

if isfield(geometry, 'tetgen_options')
    assert(isstring(geometry.tetgen_options) || ischar(geometry.tetgen_options));
    geometry = rmfields(geometry, {'refinement'});
elseif ~isfield(geometry, 'refinement')
    geometry.refinement = -1;
end
end

function geometry = check_ecs(geometry)
    if isfield(geometry, "ecs_shape")
        assert(ismember(geometry.ecs_shape, ["no_ecs", "box", "convex_hull", "tight_wrap"]));
    else
        geometry.ecs_shape = "no_ecs";
    end

    if ismember(geometry.ecs_shape, ["box", "convex_hull", "tight_wrap"])
        assert(geometry.ecs_ratio > 0);
    else
        geometry = rmfields(geometry, {'ecs_ratio'});
    end
end

function geometry = check_in(geometry)
    if isfield(geometry, "include_in") && geometry.include_in
        assert(0 < geometry.in_ratio && geometry.in_ratio < 1);
    else
        geometry.include_in = false;
    end
    if ~geometry.include_in
        geometry = rmfields(geometry, {'in_ratio'});
    end
end

function geometry = check_deformation(geometry)
    if ~isfield(geometry, 'deformation')
        geometry.deformation=[0, 0];
    end
    assert(numel(geometry.deformation) == 2);
end
