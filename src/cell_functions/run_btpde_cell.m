function [btpde_cell,btpde_soma,btpde_neurites] = run_btpde_cell(femesh_cell, setup,save_magnetization,femesh_soma,femesh_neurites,include_cell)
%RUN_BTPDE_CELL Compute the solution to the BTPDE using Matrix Formalism.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP) solves the BTPDE and returns results.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP) saves the results of each iteration at
%   "saved_simul/cell/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration. The directory
%   saved_simul can be changed with setup.saved_simul_loc;
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVE_MAGNETIZATION) also omits saving
%   or loading the magnetization field if SAVE_MAGNETIZATION is set to FALSE.
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVE_MAGNETIZATION,FEMESH_SOMA) solves, saves and returns the results for the full cell 
%   and for the simulation with only the soma. The soma signal is saved in
%  "saved_simul/soma/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT",
%  and magnetization is saved if SAVE_MAGNETIZATION is set to TRUE. The directory
%   saved_simul can be changed with setup.saved_simul_loc;
%
%   RUN_BTPDE_CELL(FEMESH, SETUP, SAVE_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES) solves, saves and returns the results for the full cell 
%   and for the simulation with only the soma and for the simulation with only each neurite branch. The ith neurite signal is saved in
%  "saved_simul/neurite_%i/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT",
%  and magnetization is saved if SAVE_MAGNETIZATION is set to TRUE. The directory
%   saved_simul can be changed with setup.saved_simul_loc;
%
%   RUN_BTPDE_CELL(FEMESH, SETUP,
%   SAVE_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES,INCLUDE_CELL) if
%   INCLUDE_CELL is set to FALSE then only the soma and neurite signals are
%   computed and saved.
%
%
%   femesh_cell: struct
%   setup: struct
%   save_magnetization (optional): logical. Defaults to true.
%   femesh_soma (optional): struct 
%   femesh_neurites (optional): struct
%   include_cell (optional): logical. Defaults to true.

include_cell = nargin < 6 || include_cell;
include_soma = nargin >= 4 && isstruct(femesh_soma);
include_neurites = nargin >= 5 && (isstruct(femesh_neurites) || iscell(femesh_neurites));
save_magnetization = (nargin <2) || save_magnetization;

if isfield(setup,'saved_simul_loc')
    savepath_root= create_savepath(setup, "btpde",setup.saved_simul_loc);
else
    savepath_root= create_savepath(setup, "btpde");
end

fprintf("Simulations to be stored in:\n%s\n",savepath_root);

if include_cell
    savepath_cell = sprintf("%s/cell",savepath_root);
    btpde_cell = solve_btpde(femesh_cell, setup,savepath_cell,save_magnetization);
else
    btpde_cell = "Not assigned";
end

if include_soma
    savepath_soma = sprintf("%s/cell",savepath_root);
    btpde_soma = solve_btpde(femesh_soma, setup,savepath_soma,save_magnetization);
else
    btpde_soma = "Not assigned";
end

if include_neurites
    nneurites=length(femesh_neurites);
    btpde_neurites = cell(nneurites,1);
    for i = 1:nneurites
        savepath_neurite= sprintf("%s/neurite_%d",savepath_root,i);
        btpde_neurites{i} = solve_btpde(femesh_neurites{i}, setup,savepath_neurite,save_magnetization);
    end
else
    btpde_neurites = "Not assigned";
end
