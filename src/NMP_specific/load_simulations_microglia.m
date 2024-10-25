function [results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_microglia(mesh,setup_file,tetgen_options,swc_file,soma_file)
%LOAD_SIMULATIONS_MICROGLIA loads experiments from a microglia cell.
%   
%   LOAD_SIMULATIONS_MICROGLIA(MESH,SETUP_FILE) loads the signal from the
%   setup file for the mesh file only, with tetgen options specified in the
%   setup file.
% 
%   LOAD_SIMULATIONS_MICROGLIA(MESH,SETUP_FILE,TETGEN_OPTIONS) loads the signal from the
%   setup file for the mesh file only, with tetgen options specified in the
%   setup file.
% 
%   LOAD_SIMULATIONS_MICROGLIA(MESH,SETUP_FILE,TETGEN_OPTIONS,SWC_FILE,SOMA_FILE) loads the signal from the
%   setup file for the soma and process of the microglia cell only, with tetgen options specified in the
%   setup file.
% 
%   mesh : (str) path to mesh
%   setup_file: (str) name of setup file.
%   tetgen_options : (str) (optional) tetgen_options for mesh. Defaults to value from
%                   setup_file.
%   swc_file : (str)(optional) path to swc file.
%   soma_file : (str) (optional) path to soma file.

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

segment_cell = nargin == 5;
if segment_cell
    setup.cell.swc = swc_file; setup.cell.soma = soma_file;
end


[setup, femesh_cell, ~, ~,femesh_soma,femesh_neurites]  = prepare_simulation(setup);
savepath_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);



if segment_cell
    max_volume = 1;
    disp("Loading simulations for soma and neurites/processes only");
    nneurites = length(femesh_neurites);
    volumes = zeros(nneurites,1);
    for i =1:nneurites
    volumes(i)= femesh_neurites{i}.total_volume;
    end
    [sorted_volumes,idx_neurites] =sort(volumes);
    idx_neurites = idx_neurites(sorted_volumes > max_volume);
    femesh_neurites = femesh_neurites(idx_neurites);
    
    load_magnetization = false; 
    savepath_soma = sprintf("%s/soma",savepath_root);
    if isfield(setup,'mf')
        mf_soma = load_mf(setup,savepath_soma,load_magnetization);
        mf_neurites = cell(nneurites,1);
        for ib = 1:nneurites
            mf_neurites{ib} =  load_mf(setup,sprintf("%s/neurite_%d",savepath_root,ib),load_magnetization);
        end
        mf_neurites = mf_neurites(idx_neurites);
        results.mf_soma = mf_soma; results.mf_neurites = mf_neurites;
    end

    if isfield(setup,'btpde')
        btpde_soma = load_btpde(setup,savepath_soma,load_magnetization);
        btpde_neurites = cell(nneurites,1);
        for ib = 1:nneurites
            btpde_neurites{ib} =  load_btpde(setup,sprintf("%s/neurite_%d",savepath_root,ib),load_magnetization);
        end
        btpde_neurites = btpde_neurites(idx_neurites);

        results.btpde_soma = btpde_soma; results.btpde_neurites = btpde_neurites;
    end


else
    disp("Loading simulations for cell only");
    savepath_cell = sprintf("%s/cell",savepath_root);
    load_magnetization=false;
    if isfield(setup,'mf')
        mf_cell= load_mf(setup,savepath_cell,load_magnetization);
       
        results.mf_cell = mf_cell;
    end

    if isfield(setup,'btpde')
        btpde_cell= load_btpde(setup,savepath_cell,load_magnetization);
       
        results.btpde_cell = btpde_cell; 
    end
end

results.setup = setup; 
toc