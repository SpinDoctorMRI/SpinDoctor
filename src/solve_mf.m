function results = solve_mf(femesh, setup, lap_eig, savepath, save_magnetization)
%SOLVE_MF Compute the solution to the BTPDE using Matrix Formalism.
%
%   SOLVE_MF(FEMESH, SETUP, LAP_EIG) solves the BTPDE using
%   Matrix Formalism and returns results.
%
%   SOLVE_MF(FEMESH, SETUP, LAP_EIG, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   SOLVE_MF(FEMESH, SETUP, LAP_EIG, SAVEPATH, SAVE_MAGNETIZATION) also omits saving
%   or loading the magnetization field if SAVE_MAGNETIZATION is set to FALSE.
%
%   femesh: struct
%   setup: struct
%   lap_eig: struct with fields
%       values: double(neig, 1)
%       funcs: double(npoint, neig)
%   savepath (optional): path string
%   save_magnetization (optinal): logical. Defaults to true.
%
%   results: struct with fields
%       magnetization: {ncompartment x namplitude x nsequence x
%                       ndirection}[npoint x 1]
%           Magnetization field at final timestep
%       signal: [ncompartment x namplitude x nsequence x ndirection]
%           Compartmentwise total magnetization at final timestep
%       signal_allcmpts: [namplitude x nsequence x ndirection]
%           Total magnetization at final timestep
%       itertimes: [namplitude x nsequence x ndirection]
%           Computational time for each iteration
%       totaltime: [1 x 1]
%           Total computational time, including matrix assembly

if numel(setup.mf.gpu) > 1
    results = solve_mf_gpus(femesh, setup, lap_eig, savepath, save_magnetization);
elseif setup.mf.gpu
    results = solve_mf_gpu(femesh, setup, lap_eig, savepath, save_magnetization);
else
    results = solve_mf_cpu(femesh, setup, lap_eig, savepath, save_magnetization);
end
