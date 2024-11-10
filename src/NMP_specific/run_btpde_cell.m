function [btpde_cell,btpde_soma,btpde_neurites] = run_btpde_cell(femesh_cell, setup, savepath_root,save_magnetization,femesh_soma,femesh_neurites,include_cell)
%RUN_BTPDE_CELL Compute the solution to the BTPDE using Matrix Formalism.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP) solves the BTPDE and returns results.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT) saves the results of each iteration at
%   "<SAVEPATH_ROOT>/cell/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT, SAVE_MAGNETIZATION) also omits saving
%   or loading the magnetization field if SAVE_MAGNETIZATION is set to FALSE.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT, SAVE_MAGNETIZATION,FEMESH_SOMA) solves, saves and returns the results for the full cell 
%   and for the simulation with only the soma. The soma signal is saved in
%  "<SAVEPATH_ROOT>/soma/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT",
%  and magnetization is saved if SAVE_MAGNETIZATION is set to TRUE.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT, SAVE_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES) solves, saves and returns the results for the full cell 
%   and for the simulation with only the soma and for the simulation with only each neurite branch. The ith neurite signal is saved in
%  "<SAVEPATH_ROOT>/neurite_%i/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT",
%  and magnetization is saved if SAVE_MAGNETIZATION is set to TRUE.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVEPATH_ROOT,
%   SAVE_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES,INCLUDE_CELL) if
%   INCLUDE_CELL is set to FALSE then only the soma and neurite signals are
%   computed and saved.
%
%
%   femesh_cell: struct
%   setup: struct
%   savepath (optional): string
%   save_magnetization (optional): logical. Defaults to true.
%   femesh_soma (optional): struct 
%   femesh_neurites (optional): struct
%   include_cell (optional): logical. Defaults to true.

include_cell = nargin < 7 || include_cell;
include_soma = nargin >= 5 && isstruct(femesh_soma);
include_neurites = nargin >= 6 && isstruct(femesh_neurites);



if nargin == 2
    btpde_cell = solve_btpde(femesh_cell, setup);
elseif nargin == 3
    save_path_cell = sprintf("%s/cell",savepath_root);
    btpde_cell = solve_btpde(femesh_cell, setup,save_path_cell);
elseif nargin >= 4 && nargin <= 6
    save_path_cell = sprintf("%s/cell",savepath_root);
    btpde_cell = solve_btpde(femesh_cell, setup,save_path_cell,save_magnetization);
elseif nargin == 7 && include_cell
    save_path_cell = sprintf("%s/cell",savepath_root);
    btpde_cell = solve_btpde(femesh_cell, setup,save_path_cell,save_magnetization);
else 
    btpde_cell = "Not assigned";
end

if include_soma
    save_path_soma = sprintf("%s/soma",savepath_root);
    btpde_soma = solve_btpde(femesh_soma, setup,save_path_soma,save_magnetization);
else
    btpde_soma = "Not assigned";
end

if include_neurites
    ndendrites=length(femesh_neurites);
    btpde_neurites = cell(ndendrites,1);
    for i=1:ndendrites 
        save_path_neurite= sprintf("%s/neurite_%d",savepath_root,i);
        btpde_neurites{i} = solve_btpde(femesh_neurites{i}, setup,save_path_neurite,save_magnetization);
    end
else
    btpde_neurites = "Not assigned";
end