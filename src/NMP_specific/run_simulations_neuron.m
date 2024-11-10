function [results,femesh_cell,femesh_soma,femesh_neurites]= run_simulations_neuron(mesh,setup_file,tetgen_options,swc_file)
%RUN_SIMULATIONS_NEURON computes experiments from a neuron.
%   
%   RUN_SIMULATIONS_NEURON(MESH,SETUP_FILE) computes the signal from the
%   setup file for the mesh file only, with tetgen options specified in the
%   setup file.
% 
%   RUN_SIMULATIONS_NEURON(MESH,SETUP_FILE,TETGEN_OPTIONS) computes the signal from the
%   setup file for the mesh file only, with tetgen options specified in the
%   setup file.
% 
%   RUN_SIMULATIONS_NEURON(MESH,SETUP_FILE,TETGEN_OPTIONS,SWC_FILE) computes the signal from the
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
[~,cellname,~] = fileparts(setup.name);

segment_cell = nargin == 4;
if segment_cell
    setup.cell.swc = swc_file;
end

[setup, femesh_cell, ~, ~,femesh_soma,femesh_neurites]  = prepare_simulation(setup);
if isfield(setup,'saved_simul_loc')
savepath_root=create_savepath(setup,"mf",setup.saved_simul_loc);
else
savepath_root=create_savepath(setup,"mf");
end

save_magnetization = false; 

if segment_cell
    disp("Running simulations for soma and neurites/processes only");
    include_cell = false;
    if isfield(setup,'mf')
        [mf_cell,mf_soma,mf_neurites,~,~,~] = run_mf_cell(...
            femesh_cell, setup, savepath_root,save_magnetization,femesh_soma,femesh_neurites,include_cell);
        results.mf_cell = mf_cell; results.mf_soma = mf_soma;results.mf_neurites = mf_neurites;
    end
    
    if isfield(setup,'btpde')
       [btpde_cell,btpde_soma,btpde_neurites] = run_btpde_cell(...
           femesh_cell, setup, savepath_root,save_magnetization,femesh_soma,femesh_neurites,include_cell);
        results.btpde_cell = btpde_cell; results.btpde_soma = btpde_soma;results.btpde_neurites = btpde_neurites;
    end
else
    disp("Running simulations for cell only");
    if isfield(setup,'mf')
        [mf_cell,~,~,~,~,~] = run_mf_cell(...
            femesh_cell, setup, savepath_root,save_magnetization);
        results.mf_cell = mf_cell; 
    end
    if isfield(setup,'btpde')
       [btpde_cell,~,~] = run_btpde_cell(...
           femesh_cell, setup, savepath_root,save_magnetization);
        results.btpde_cell = btpde_cell; 
    end
end

results.setup = setup; 
toc