function [results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_neuron(mesh,setup_file,tetgen_options,swc_file)
%LOAD_SIMULATIONS_NEURON loads experiments from a neuron.
%   
%   LOAD_SIMULATIONS_NEURON(MESH,SETUP_FILE) loads the signal from the
%   setup file for the mesh file only, with tetgen options specified in the
%   setup file.
% 
%   LOAD_SIMULATIONS_NEURON(MESH,SETUP_FILE,TETGEN_OPTIONS) loads the signal from the
%   setup file for the mesh file only, with tetgen options specified in the
%   setup file.
% 
%   LOAD_SIMULATIONS_NEURON(MESH,SETUP_FILE,TETGEN_OPTIONS,SWC_FILE) loads the signal from the
%   setup file for the soma and process of the microglia cell only, with tetgen options specified in the
%   setup file.
% 
%   mesh : (str) path to mesh
%   setup_file: (str) name of setup file.
%   tetgen_options : (str) (optional) tetgen_options for mesh. Defaults to value from
%                   setup_file.
%   swc_file : (str)(optional) path to swc file.

tic
addpath(genpath('setups'));
addpath(genpath('src'));

% Launch setup here
run(sprintf("%s.m",setup_file));
fprintf("Running %s.m\n",setup_file)


if nargin >=3
    setup.geometry.tetgen_options=string(tetgen_options);
else
    fprintf('Applying default Tetgen parameters %s\n',setup.geometry.tetgen_options);
end

setup.name=string(mesh);
% [~,cellname,~] = fileparts(setup.name);

segment_cell = nargin == 4;

if segment_cell
    setup.cell.swc = swc_file;
end


[setup, femesh_cell, ~, ~,femesh_soma,femesh_neurites]  = prepare_simulation(setup);

if segment_cell
    disp("Loading simulations for soma and neurites/processes only");
    nneurites = length(femesh_neurites); 
    include_soma = true; include_cell = false;
else    
    disp("Loading simulations for cell only");
    nneurites = 0; include_soma = false; include_cell = true;
end

if isfield(setup,'mf')
    [mf_cell,mf_soma,mf_neurites] = load_mf_cell(setup,false,nneurites,include_soma,include_cell);
    results.mf_cell = mf_cell; results.mf_soma = mf_soma;results.mf_neurites = mf_neurites;
end


if isfield(setup,'btpde')
    [btpde_cell,btpde_soma,btpde_neurites] = load_btpde_cell(setup,false,nneurites,include_soma,include_cell);
    results.btpde_cell = btpde_cell; results.btpde_soma = btpde_soma;results.btpde_neurites = btpde_neurites;
end
results.setup = setup; 
toc