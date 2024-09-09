function [femesh, surfaces, cells] = create_geometry(setup)
%CREATE_GEOMETRY Create cells, surfaces and finite element mesh.
%   This function does the following:
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
%       nboundary: [1 x 1]
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
cell_shape = setup.geometry.cell_shape;
ncompartment = setup.ncompartment;
nboundary = setup.nboundary;

% Check that folder exists
[filepath, name, ext] = fileparts(filename);
if ~isfolder(filepath) && ~isempty(filepath)
    mkdir(filepath);
end
% is public polyhedral file formats (STL, PLY, OFF)
is_ply = ismember(lower(ext), [".stl", ".ply", ".off"]);

% File name to save or load cell description
cellfilename = filename + "_cells";

% Check if cell description file is already available
if isfile(cellfilename)
    cells = read_cells(cellfilename);
elseif ismember(cell_shape, ["sphere" "cylinder"])
    % Create cells
    cells = create_cells(setup);
    % Save cell configuration
    save_cells(cells, cellfilename);
else
    % No cells
    cells = struct;
end

% Make directory for storing finite elements mesh
meshdir_name = name + replace(ext, '.', '_') + "_dir";
save_meshdir_path = fullfile(filepath, meshdir_name);
if ~isfolder(save_meshdir_path)
    mkdir(save_meshdir_path);
end
if is_ply
    filename_new = save_meshdir_path + "/" + name + ext;
    if ~isfile(filename_new)
        copyfile(filename, save_meshdir_path);
        if isfield(setup.geometry, "tetgen_options")
            call_tetgen(filename_new, setup.geometry.tetgen_options);
        else    
            call_tetgen(filename_new);
        end
    end
    tetgen_ply = replace(filename_new, ext, '.1');
end

% Use an existing finite elements mesh. 
% The name of the finite elements mesh is stored in the string fname_tetgen_femesh
refinement_str = "";
if isfield(setup.geometry, "refinement")
    refinement_str = sprintf("_refine%g", setup.geometry.refinement);
end
if isfield(setup.geometry, "tetgen_options")
    refinement_str = sprintf("_tet%s", setup.geometry.tetgen_options);
end
ecs_str = sprintf("_%s", setup.geometry.ecs_shape);
if isfield(setup.geometry, 'ecs_ratio') && setup.geometry.ecs_shape ~= "no_ecs"
    ecs_str = ecs_str + sprintf("%g", setup.geometry.ecs_ratio);
end

fname_tetgen = save_meshdir_path + "/" + name + ecs_str + refinement_str + "_mesh";

% Read or create surface triangulation
if isfile(fname_tetgen + ".node") && isfile(fname_tetgen + ".poly")
    surfaces = read_surfaces(fname_tetgen);
else
    switch cell_shape
        case "sphere"
            % Create surface geometry of spheres
            surfaces = create_surfaces_sphere(cells, setup);
        case "cylinder"
            % Create surface geometry of cylinders
            surfaces = create_surfaces_cylinder(cells, setup);
        case "neuron"
            if is_ply
                surfaces = create_surfaces_neuron(tetgen_ply, setup);
            else
                surfaces = create_surfaces_neuron(filename, setup);
            end
    end
    disp(fname_tetgen)
    save_surfaces(fname_tetgen, surfaces);
end

% Add ".1" suffix to output file name, since this is what Tetgen does
fname_tetgen_femesh = fname_tetgen + ".1";

if isfield(setup.geometry, "tetgen_options")
    tetgen_params = {setup.geometry.tetgen_options};
elseif isfield(setup.geometry, "refinement")
    tetgen_params = {setup.geometry.refinement};
else
    tetgen_params = {};
end

if ~isfile(fname_tetgen_femesh + ".node")
    call_tetgen(fname_tetgen + ".poly", tetgen_params{:});
end

try
    % Read global mesh from Tetgen output
    femesh_all = read_tetgen(fname_tetgen_femesh);

    % Check that at correct number of compartments and boundaries has been found
    compartments = unique(femesh_all.elementmarkers);
    boundaries = unique(femesh_all.facetmarkers);
    solution = "use smaller refinement or change surface triangulation.";
    assert(ncompartment == length(compartments), "Incorrect number of compartments, " + solution);
    assert(nboundary == length(boundaries), "Incorrect number of boundaries, " + solution);
catch ME
    delete(fname_tetgen+"*")
    rethrow(ME)
end

% Deform domain
if any(setup.geometry.deformation)
    fprintf("Deforming domain with bend %g and twist %g\n", setup.geometry.deformation);
    femesh_all.points = deform_domain(femesh_all.points, setup.geometry.deformation);
end

% Split mesh into compartments
femesh = split_mesh(femesh_all);

% Get volume and surface area quantities from mesh
[volumes, areas] = get_vol_sa(femesh);
femesh.volumes = volumes;
femesh.total_volume = sum(volumes, 'all');
femesh.areas = areas;
femesh.total_area = sum(areas, 'all');

surfaces.areas = areas;
surfaces.total_area = sum(areas, 'all');
