function [femesh, surfaces] = create_femesh(cells, params_cells, params_domain)
%CREATE_FEMESH Create finite element mesh.
%   This function first creates and saves or load a surface geometry, calls
%   Tetgen, deforms the domain, and splits mesh into compartments.
%
%   cells: struct
%   params_cells: struct
%
%   femesh: struct with fields
%       ncompartment: [1 x 1]
%     	nboundary: [1 x 1]
%       points: {1 x ncompartment}[3 x npoint(icmpt)]
%       facets: {ncompartment x nboundary}[3 x nfacet(icmpt, ibdry)]
%       elements: {1 x ncompartment}[4 x nelement(icmpt)]
%   	point_map: {1 x ncompartment}[1 x npoint(icmpt)]
%   surfaces: struct with fields
%       points: [3 x npoint]
%       facets: [3 x nfacet]
%       facetmarkers: [1 x nfacet]
%       regions: [3 x ncompartment]


% Extract parameters
filename = params_cells.filename;
shape = params_cells.shape;
ncompartment = length(params_domain.compartments);
nboundary = length(params_domain.boundaries);

% Make directory for storing finite elements mesh
is_stl =  endsWith(filename, ".stl");
if is_stl
    assert(params_cells.ecs_shape == "no_ecs");
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
if isfield(params_cells, "refinement")
    refinement_str = sprintf("_refinement%g", params_cells.refinement);
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
            surfaces = create_surfaces_sphere(cells, params_cells);
        case "cylinder"
            % Create surface geometry of cylinders
            surfaces = create_surfaces_cylinder(cells, params_cells);
        case "neuron"
            if is_stl
                surfaces = struct;
            else
                surfaces = create_surfaces_neuron(filename, params_cells);
            end
    end

    if ~is_stl
        % plot_surface_triangulation(surfaces);
        save_surfaces(fname_tetgen, surfaces);
    end
end

% Add ".1" suffix to output file name, since this is what Tetgen does
fname_tetgen_femesh = fname_tetgen + ".1";

if isfield(params_cells, "refinement")
    tetgen_params = {params_cells.refinement};
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
solution = "use smaller hmax or change surf triangulation.";
assert(ncompartment == length(compartments), "Incorrect number of compartments, " + solution);
assert(nboundary == length(boundaries), "Incorrect number of boundaries, " + solution);

% Deform domain
if any(params_cells.deformation)
    fprintf("Deforming domain with bend %g and twist %g\n", params_cells.deformation);
    femesh_all.points = deform_domain(femesh_all.points, params_cells.deformation);
end

% Split mesh into compartments
femesh = split_mesh(femesh_all);
