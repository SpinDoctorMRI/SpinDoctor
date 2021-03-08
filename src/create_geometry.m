function [femesh, surfaces, cells] = create_geometry(setup)
%CREATE_GEOMETRY Create cells, surfaces and finite element mesh.
%   This function does the following:
%   - Check geometry setup consistency
%   - Create or load cell configuration
%   - Create or load surface triangulation
%   - Call TetGen
%   - Deform domain
%   - Split mesh into compartments
%   For custom geometries with more than one compartment, call `split_mesh`
%   directly instead. This requires facet and element labels.
%
%   setup: struct
%
%   femesh: struct with fields
%       ncompartment: [1 x 1]
%         nboundary: [1 x 1]
%       points: {1 x ncompartment}[3 x npoint(icmpt)]
%       facets: {ncompartment x nboundary}[3 x nfacet(icmpt, ibdry)]
%       elements: {1 x ncompartment}[4 x nelement(icmpt)]
%       point_map: {1 x ncompartment}[1 x npoint(icmpt)]
%   surfaces: struct with fields
%       points: [3 x npoint]
%       facets: [3 x nfacet]
%       facetmarkers: [1 x nfacet]
%       regions: [3 x ncompartment]
%   cells: struct with fields
%       centers: [d x ncell], where d is 2 (cylinders) or 3 (spheres)
%       radii: [1 x ncell]


% Extract parameters
filename = setup.name;
shape = setup.geometry.cell_shape;
ncompartment = length(setup.pde.compartments);
nboundary = length(setup.pde.boundaries);


% Check correct input format
assert(ismember(setup.geometry.cell_shape, ["sphere" "cylinder" "neuron"]))
if setup.geometry.cell_shape == "neuron"
    % We do not deform the finite elements mesh if in the Neuron module. We do
    % deform the finite elements mesh if in SpinDoctor White Matter mode.
    if isfield(setup.geometry, "deformation") && any(setup.geometry.deformation ~= 0)
        error("Deformation is not available for neurons. Set deformation to [0; 0].")
    end
    if setup.geometry.ncell ~= 1
        error("Neuron cell type is only available for ncell=1.")
    end
end
assert(~(setup.geometry.ecs_shape == "no_ecs" && setup.geometry.ncell > 1), ...
    "Geometry must include ECS if more than one cell");
assert(ismember(setup.geometry.ecs_shape, ["no_ecs" "box" "convex_hull" "tight_wrap"]))


% Check that folder exists
parts = split(filename, "/");
if length(parts) >= 2
    folder = join(parts(1:end-1), "/");
    if ~isfolder(folder)
        mkdir(folder);
    end
end


% File name to save or load cell description
cellfilename = filename + "_cells";

% Check if cell description file is already available
if isfile(cellfilename)
    cells = read_cells(cellfilename);
elseif ismember(setup.geometry.cell_shape, ["sphere" "cylinder"])
    % Create cells
    cells = create_cells(setup);
    
    % Save cell configuration
    save_cells(cells, cellfilename);
else
    % No cells
    cells = struct;
end


% Make directory for storing finite elements mesh
is_stl = endsWith(filename, ".stl");
if is_stl
    % ECS is currently only available for surface meshes
    assert(setup.geometry.ecs_shape == "no_ecs");
    tmp = split(filename, ".stl");
    filename = tmp(1);
end
save_meshdir_path = filename + "_dir";
if ~isfolder(save_meshdir_path)
    mkdir(save_meshdir_path);
end


% Use an existing finite elements mesh or create a new finite
% elements mesh. The name of the finite elements mesh is stored in the string
% fname_tetgen_femesh
refinement_str = "";
if isfield(setup.geometry, "refinement")
    refinement_str = sprintf("_refinement%g", setup.geometry.refinement);
end
tmp = split(filename, "/");
fname_tetgen = save_meshdir_path + "/" + tmp(end) + refinement_str + "_mesh";

% Read or create surface triangulation
if isfile(fname_tetgen + ".node") && isfile(fname_tetgen + ".poly")
    surfaces = read_surfaces(fname_tetgen);
else
    switch shape
        case "sphere"
            % Create surface geometry of spheres
            surfaces = create_surfaces_sphere(cells, setup);
        case "cylinder"
            % Create surface geometry of cylinders
            surfaces = create_surfaces_cylinder(cells, setup);
        case "neuron"
            if is_stl
                surfaces = struct;
            else
                surfaces = create_surfaces_neuron(filename, setup);
            end
    end

    if ~is_stl
        % plot_surface_triangulation(surfaces);
        save_surfaces(fname_tetgen, surfaces);
    end
end

% Add ".1" suffix to output file name, since this is what Tetgen does
fname_tetgen_femesh = fname_tetgen + ".1";

if isfield(setup.geometry, "refinement")
    tetgen_params = {setup.geometry.refinement};
else
    tetgen_params = {};
end

if isfile(fname_tetgen_femesh + ".node")
elseif is_stl
    call_tetgen(filename + ".stl", tetgen_params{:});
    fname_tetgen_femesh = filename + ".1";
else
    call_tetgen(fname_tetgen + ".poly", tetgen_params{:});
end

% Read global mesh from Tetgen output
femesh_all = read_tetgen(fname_tetgen_femesh);

% Check that at correct number of compartments and boundaries has been found
compartments = unique(femesh_all.elementmarkers);
boundaries = unique(femesh_all.facetmarkers);
solution = "use smaller refinement or change surface triangulation.";
assert(ncompartment == length(compartments), "Incorrect number of compartments, " + solution);
assert(nboundary == length(boundaries), "Incorrect number of boundaries, " + solution);

% Deform domain
if any(setup.geometry.deformation)
    fprintf("Deforming domain with bend %g and twist %g\n", setup.geometry.deformation);
    femesh_all.points = deform_domain(femesh_all.points, setup.geometry.deformation);
end

% Split mesh into compartments
femesh = split_mesh(femesh_all);
