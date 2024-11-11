function [btpde_cell,btpde_soma,btpde_neurites] = load_btpde_cell(setup, savepath_root,load_magnetization,nneurites,include_soma,include_cell)
%LOAD_BTPDE_CELL Compute the solution to the BTPDE using Matrix Formalism.
%
%   LOAD_BTPDE_CELL(FEMESH, SETUP) solves the BTPDE and returns results.
%
%   LOAD_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT) loads the results of each iteration from
%   "<SAVEPATH_ROOT>/cell/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%
%   LOAD_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT, LOAD_MAGNETIZATION) also omits loading the magnetization field if LOAD_MAGNETIZATION is set to FALSE.
%
%   LOAD_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT, LOAD_MAGNETIZATION,FEMESH_SOMA) eturns the results for the full cell 
%   and for the simulation with only the soma. The soma signal is taken from
%  "<SAVEPATH_ROOT>/soma/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT",
%  and magnetization is loaded if LOAD_MAGNETIZATION is set to TRUE.
%
%   LOAD_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT, LOAD_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES)  returns the results for the full cell 
%   and for the simulation with only the soma and for the simulation with only each neurite branch. The ith neurite signal is from
%  "<SAVEPATH_ROOT>/neurite_%i/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT",
%  and magnetization is loaded if LOAD_MAGNETIZATION is set to TRUE.
%
%   LOAD_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT,
%   LOAD_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES,INCLUDE_CELL) if
%   INCLUDE_CELL is set to FALSE then only the soma and neurite signals are
%   loaded.
%
%
%   femesh_cell: struct
%   setup: struct
%   savepath (optional): string
%   load_magnetization (optional): logical. Defaults to true.
%   femesh_soma (optional): struct 
%   femesh_neurites (optional): struct
%   include_cell (optional): logical. Defaults to true.

if include_cell
    save_path_cell = sprintf("%s/cell",savepath_root);
    btpde_cell = load_btpde( setup,save_path_cell);
elseif nargin >= 4 && nargin <= 6
    save_path_cell = sprintf("%s/cell",savepath_root);
    btpde_cell = load_btpde(setup,save_path_cell,load_magnetization);
elseif nargin == 7 && include_cell
    save_path_cell = sprintf("%s/cell",savepath_root);
    btpde_cell = load_btpde( setup,save_path_cell,load_magnetization);
else 
    btpde_cell = "Not assigned";
end

if nargin >= 5
    save_path_soma = sprintf("%s/soma",savepath_root);
    btpde_soma = load_btpde(setup,save_path_soma,load_magnetization);
else
    btpde_soma = "Not assigned";
end

if nargin >= 6
    ndendrites=length(femesh_neurites);
    btpde_neurites = cell(ndendrites,1);
    for i=1:ndendrites 
        save_path_neurite= sprintf("%s/neurite_%d",savepath_root,i);
        btpde_neurites{i} = load_btpde( setup,save_path_neurite,load_magnetization);
    end
else
    btpde_neurites = "Not assigned";
end