function cells = create_cells(params_cells)
%CREATE_CELLS Create cells in their canonical configuration.
%
%   params_cells: struct
%
%   cells: struct



assert(any(params_cells.shape == ["sphere", "cylinder", "neuron"]))
if params_cells.shape == "neuron"
    % We do not deform the finite elements mesh if in the Neuron module. We do
    % deform the finite elements mesh if in SpinDoctor White Matter mode.
    if isfield(params_cells, "deformation") && any(params_cells.deformation ~= 0)
        error("Deformation is not available for neurons. Set deformation to [0; 0].")
    end
    if params_cells.ncell ~= 1
        error("Neuron cell type is only available for ncell=1.")
    end
end

assert(~(params_cells.ecs_shape == "no_ecs" && params_cells.ncell > 1), ...
    "Geometry must include ECS if more than one cell");

assert(any(params_cells.ecs_shape == ["no_ecs", "box", "convex_hull", "tight_wrap"]))

% File name to save or load cell description
cellfilename = params_cells.filename + "_cells";

% Return if cell description file is already available
if isfile(cellfilename)
    cells = read_cells_inputfile(cellfilename);
    return
end

% Create cell description files
switch params_cells.shape
    case {"sphere", "cylinder"}
        cells = create_cells_inputfile(params_cells, cellfilename);
    case "neuron"
        cells = struct;
end
