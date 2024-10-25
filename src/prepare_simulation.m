function [setup, femesh, surfaces, cells,femesh_soma,femesh_neurites] = prepare_simulation(setup)
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

% Segment the cell, if options are enabled.
if isfield(setup,'cell') && isfield(setup.cell,'swc')
    [mesh_path,cellname,~] = fileparts(setup.name);
    tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
    if isfile(setup.cell.swc)
        swc_file = setup.cell.swc;
        if isfield(setup.cell,'soma')
            soma_file = setup.cell.soma;
            if not(isfile(soma_file))
                error("File %s was used to segment the mesh but no such file exists.\n",soma_file);
            end
           [femesh,femesh_soma,femesh_neurites] = segment_femesh(femesh,swc_file,tetgen_path,soma_file); 
        else
           [femesh,femesh_soma,femesh_neurites] = segment_femesh(femesh,swc_file,tetgen_path); 
        end
    else
        disp('error')
        error("File %s was used to segment the mesh but no such file exists.\n",setup.cell.swc);
    end
else
    femesh_soma = "Not assigned"; femesh_neurites = "Not assigned";
end
