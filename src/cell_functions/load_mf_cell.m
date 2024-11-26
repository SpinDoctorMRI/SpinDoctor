function [mf_cell,mf_soma,mf_neurites] = load_mf_cell(setup,load_magnetization,nneurites,include_soma,include_cell)
%LOAD_MF_CELL Compute the solution to the BTPDE using Matrix Formalism.
%
%
%   LOAD_MF_CELL(SETUP) loads the results of for the cell. The results are taken from
%   "saved_simul/cell/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   The directory saved_simul can be changed with setup.saved_simul_loc;

%   LOAD_MF_CELL(SETUP, LOAD_MAGNETIZATION) also omits 
%   loading the magnetization field if LOAD_MAGNETIZATION is set to FALSE.
%
%   LOAD_MF_CELL(SETUP, LOAD_MAGNETIZATION,NNEURITES) 
%   loads the results for the neurites from 1 to
%   NNEURITES.
% 
%   LOAD_MF_CELL(SETUP, SAVEPATH_ROOT,
%   LOAD_MAGNETIZATION,NNEURITES,INCLUDE_SOMA) loads the results for the 
%   soma if INCLUDE_SOMA is true and for the neurites from 1 to NNEURITES.
% 
%   LOAD_MF_CELL(SETUP, SAVEPATH_ROOT,
%   LOAD_MAGNETIZATION,NNEURITES,INCLUDE_SOMA) loads the results for the 
%   cell if INCLUDE_CELL is true,
%   the soma if INCLUDE_SOMA is true 
%   and for the neurites from 1 to NNEURITES.

%   setup: struct
%   load_magnetization (optional): logical. Defaults to false.
%   nneurites (optional): int. Defaults to 0.
%   include_soma (optional): logical. Defaults to false
%   include_cell (optional): logical. Defaults to false.


load_magnetization = nargin >= 2 && load_magnetization; 
include_cell = nargin<=2 || (nargin ==5 && include_cell);
include_soma = nargin >=4 && include_soma;
if nargin < 3
    nneurites = 0;
end

savepath_root= create_savepath(setup, "mf");

fprintf("Simulations to be stored in:\n%s\n",savepath_root);

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
if nneurites > 0
    for ib = 1:nneurites
        savepath_neurite = sprintf('%s/neurite_%d',savepath_root,ib);
        mf_neurites{ib} = load_mf(setup,savepath_neurite,load_magnetization);
    end
else
    mf_neurites = "Not assigned";
end


