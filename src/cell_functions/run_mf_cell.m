function [mf_cell,mf_soma,mf_neurites,lap_eig_cell,lap_eig_soma,lap_eig_neurites] = run_mf_cell(femesh_cell, setup,save_magnetization,femesh_soma,femesh_neurites,include_cell)
%RUN_MF_CELL Compute the solution to the BTPDE using Matrix Formalism.
%
%
%   RUN_MF_CELL(FEMESH, SETUP) saves the results of each iteration at
%   "saved_simul/cell/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration. The directory
%   saved_simul can be changed with setup.saved_simul_loc;
%
%   RUN_MF_CELL(FEMESH, SETUP, SAVE_MAGNETIZATION) also omits saving
%   or loading the magnetization field if SAVE_MAGNETIZATION is set to FALSE.
%
%   RUN_MF_CELL(FEMESH, SETUP, SAVE_MAGNETIZATION,FEMESH_SOMA) solves, saves and returns the results for the full cell 
%   and for the simulation with only the soma. The soma signal is saved in
%  "saved_simul/soma/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT",
%  and magnetization is saved if SAVE_MAGNETIZATION is set to TRUE.The directory
%   saved_simul can be changed with setup.saved_simul_loc;
%
%   RUN_MF_CELL(FEMESH, SETUP, SAVE_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES) solves, saves and returns the results for the full cell 
%   and for the simulation with only the soma and for the simulation with only each neurite branch. The ith neurite signal is saved in
%  "saved_simul/neurite_%i/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT",
%  and magnetization is saved if SAVE_MAGNETIZATION is set to TRUE. The directory
%   saved_simul can be changed with setup.saved_simul_loc;
%
%   RUN_MF_CELL(FEMESH, SETUP,
%   SAVE_MAGNETIZATION,FEMESH_SOMA,FEMESH_NEURITES,INCLUDE_CELL) if
%   INCLUDE_CELL is set to FALSE then only the soma and neurite signals are
%   computed and saved.
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

save_eig = setup.mf.save_eig;

if not(save_eig) && save_magnetization
    error("Cannot have setup.mf.save_eig = false and save magentization.");
end

savepath_root= create_savepath(setup, "mf");

fprintf("Simulations to be stored in:\n%s\n",savepath_root);



% Cell
if include_cell
    savepath_cell = sprintf("%s/cell",savepath_root);
    % FEM_save_path = make_FEM_save_path(setup.pde,setup.mf,false,savepath_cell);
    FEM_save_path = sprintf("%s/FEM_lap_eig.mat",add_mf_str_savepath(savepath_cell,setup.mf ));
    compute_eig = ~isfile(FEM_save_path);
    if compute_eig | save_magnetization
        lap_eig_cell = compute_laplace_eig(femesh_cell, setup.pde, setup.mf,savepath_cell,save_eig);  
        mf_cell = solve_mf(femesh_cell, setup, lap_eig_cell,savepath_cell,save_magnetization);
    else
        mf_cell  = solve_mf_cpu_no_lap(setup, savepath_cell);
        lap_eig_cell = "Not assigned";
    end
else
    mf_cell = "Not assigned";
    lap_eig_cell = "Not assigned";
end


% Soma
if include_soma
    savepath_soma= sprintf("%s/soma",savepath_root);
    % FEM_save_path = make_FEM_save_path(setup.pde,setup.mf,false,savepath_soma);
    FEM_save_path = sprintf("%s/FEM_lap_eig.mat",add_mf_str_savepath(savepath_soma,setup.mf ));
    compute_eig = ~isfile(FEM_save_path);
    if compute_eig | save_magnetization
        lap_eig_soma = compute_laplace_eig(femesh_soma, setup.pde, setup.mf,savepath_soma,save_eig);  
        mf_soma = solve_mf(femesh_soma, setup, lap_eig_soma,savepath_soma,save_magnetization);
    else
        mf_soma = solve_mf_cpu_no_lap(setup, savepath_soma);
        lap_eig_soma = "Not assigned";
    end
else
    mf_soma = "Not assigned";
    lap_eig_soma = "Not assigned";
end


% Neurites
if include_neurites
    nneurites = length(femesh_neurites);
    mf_neurites = cell(nneurites,1);
    lap_eig_neurites = cell(nneurites,1);
    for ib = 1:nneurites
        savepath_neurite= sprintf("%s/neurite_%d",savepath_root,ib);
        % FEM_save_path = make_FEM_save_path(setup.pde,setup.mf,false,savepath_neurite);
        FEM_save_path = sprintf("%s/FEM_lap_eig.mat",add_mf_str_savepath(savepath_neurite,setup.mf ));

        compute_eig = ~isfile(FEM_save_path);
        if compute_eig | save_magnetization
             lap_eig_neurites{ib}= compute_laplace_eig(femesh_neurites{ib}, setup.pde, setup.mf,savepath_neurite,save_eig);  
            mf_neurites{ib} = solve_mf(femesh_neurites{ib}, setup,  lap_eig_neurites{ib},savepath_neurite,save_magnetization);
        else
            mf_neurites{ib} = solve_mf_cpu_no_lap( setup, savepath_neurite);
            lap_eig_neurites = "Not assigned";
        end
    end

else
    mf_neurites = "Not assigned";
    lap_eig_neurites = "Not assigned";
end



