function [mf_cell,mf_soma,mf_neurites] = load_mf_cell(setup, savepath_root,load_magnetization,nneurites,include_soma,include_cell)
%LOAD_MF_CELL Compute the solution to the BTPDE using Matrix Formalism.
%
%
%   LOAD_MF_CELL(FEMESH, SETUP, SAVEPATH_ROOT) loads the results of each
%   iteration for the full cell and soma, the results are taken from
%   "<SAVEPATH_ROOT>/cell/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   and
%   "<SAVEPATH_ROOT>/soma/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   is loaded

%   LOAD_MF_CELL(FEMESH, SETUP, SAVEPATH, LOAD_MAGNETIZATION) also omits 
%   loading the magnetization field if LOAD_MAGNETIZATION is set to FALSE.
%
%   LOAD_MF_CELL(FEMESH, SETUP, SAVEPATH_ROOT,
%   LOAD_MAGNETIZATION,NNEURITES) loads the results of each
%   iteration for the full cell, soma, and the neurites from 1 to
%   NNEURITES.
% 
%   LOAD_MF_CELL(SETUP, SAVEPATH_ROOT,
%   LOAD_MAGNETIZATION,INCLUDE_CELL) if
%   INCLUDE_CELL is set to FALSE then only the soma and neurite signals are
%   loaded.
% 
%   setup: struct
%   savepath_root: string
%   load_magnetization (optional): logical. Defaults to false.
%   nneurites (optional): int
%   include_cell (optional): logical. Defaults to true.

load_magnetization = nargin >= 3 && load_magnetization; 

% Cell
if include_cell
    savepath_cell = sprintf("%s/cell",savepath_root);
    mf_cell = load_mf(setup,savepath_cell,load_magnetization);
else
    mf_cell = "Not assigned";
end

% Soma
if include_soma
savepath_soma = sprintf("%s/soma",savepath_root);
mf_soma = load_mf(setup,savepath_soma,load_magnetization);
else
    mf_soma = "Not assigned";
end
% Neurites
mf_neurites = cell(nneurites,1);
for ib = 1:nneurites
    savepath_neurite = sprintf('%s/neurite_%d',savepath_root,ib);
    mf_neurites{ib} = load_mf(setup,savepath_neurite,load_magnetization);
end



