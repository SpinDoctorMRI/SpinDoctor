function results = load_btpde_midpoint(setup, savepath, load_magnetization)
%LOAD_BTPDE_MIDPOINT Load the results saved by SOLVE_BTPDE_MIDPOINT.
%
%   LOAD_BTDPE_MIDPOINT(SETUP, SAVEPATH) loads the results of each iteration from
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde_midpoint/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%
%   LOAD_BTDPE_MIDPOINT(SETUP, SAVEPATH, LOAD_MAGNETIZATION) also omits loading
%   the magnetization field if LOAD_MAGNETIZATION is set to FALSE.
%
%   setup: struct
%   savepath: string
%   load_magnetization (optional): logical. Defaults to true.
%
%   results: struct with fields. Split into the experiments for constant
%   direction vector sequences (const) and those with varying direction
%   (camino) from camino files. 
%   If only const or only camino sequences are present, then this
%   additional struct layer is removed.
%   
%       camino.magnetization: {ncompartment x nsequeunce}[npoint x 1]
%          Magnetization field at final timestep
%       camino.signal: [ncompartment x nsequence]
%           Compartmentwise total magnetization at final timestep
%       camino.signal_allcmpts: [nsequence x 1]
%           Total magnetization at final timestep
%       camino.itertimes: [nsequence x 1]
%           Computational time for each iteration
%       const.magnetization: {ncompartment x namplitude x nsequence x
%                       ndirection}[npoint x 1]
%           Magnetization field at final timestep
%       const.signal: [ncompartment x namplitude x nsequence x ndirection]
%           Compartmentwise total magnetization at final timestep
%       const.signal_allcmpts: [namplitude x nsequence x ndirection]
%           Total magnetization at final timestep
%       const.itertimes: [namplitude x nsequence x ndirection]
%           Computational time for each iteration
%       totaltime: [1 x 1]
%           Total computational time, including matrix assembly


% Measure function evaluation time
starttime = tic;

% Provide default value
if nargin < nargin(@load_btpde_midpoint)
    load_magnetization = true;
end

% Extract experiment parameters
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;
theta = setup.btpde_midpoint.implicitness;
dt = setup.btpde_midpoint.timestep;

% Sizes
ncompartment = setup.ncompartment;
namplitude = setup.namplitude;
nsequence = setup.nsequence;
ndirection = setup.ndirection;

% Folder for saving
savepath = sprintf( ...
    "%s/theta%g_dt%g", ...
    savepath, theta, dt ...
);

% Initialize output arguments
const_sequences_ind = cellfun(@(x) ~isa(x,"SequenceCamino"),sequences,'UniformOutput',true);
nsequence_const = sum(const_sequences_ind);
sequences_const = sequences(const_sequences_ind);
const = struct;
const.magnetization = cell(ncompartment, namplitude, nsequence_const, ndirection);
const.signal = inf(ncompartment, namplitude, nsequence_const, ndirection);
const.signal_allcmpts = zeros(namplitude, nsequence_const, ndirection);
const.itertimes = zeros(namplitude, nsequence_const, ndirection);

nsequence_camino = sum(~const_sequences_ind);
camino = struct;
camino.magnetization = cell(ncompartment,nsequence_camino, 1);
camino.signal = inf(ncompartment, nsequence_camino);
camino.signal_allcmpts = zeros(nsequence_camino,1);
sequences_camino=sequences(~const_sequences_ind);
camino.itertimes = zeros(nsequence_camino, 1);

inds = [namplitude ndirection];
% Iterate over gradient amplitudes, sequences and directions.
% Checking for constant sequences
for iseq = 1:nsequence_const
    seq = sequences_const{iseq};
    % Load results
    filename = sprintf("%s/%s.mat", savepath, seq.string(true));    
    fprintf("Load btpde_midpoint %d/%d.\n", iseq, nsequence);
    mfile = load(filename);
    for iall = 1:prod(inds)
        [iamp, idir] = ind2sub(inds, iall);
        % Extract iteration inputs
        b = bvalues(iamp, iseq);
        ug = directions(:, idir);
        data = mfile.(gradient_fieldstring(ug, b));

        const.signal(:, iamp, iseq, idir) = data.signal;
        const.itertimes(iamp, iseq, idir) = data.itertimes;
        if load_magnetization
            const.magnetization(:, iamp, iseq, idir) = data.magnetization;
        end
    end
end

% Checking for camino sequences
for iseq = 1:nsequence_camino
    seq = sequences_camino{iseq};
    filename = sprintf("%s/%s.mat", savepath, seq.string(true));
    mfile = matfile(filename, "Writable", false);
    fprintf("Load btpde_midpoint for %s \n", seq.string);
    savedata = mfile.(seq.string);
    camino.signal(:,iseq) = savedata.signal;
    camino.itertimes(iseq) = savedata.itertimes;
    if load_magnetization
        camino.magnetization(:, iseq) = savedata.magnetization;
    end
end

% Total magnetization (sum over compartments)
camino.signal_allcmpts(:) = sum(camino.signal, 1);
const.signal_allcmpts(:) = sum(const.signal, 1);

totaltime = sum(camino.itertimes,"all") + sum(const.itertimes,"all");
results = merge_results(camino,const,nsequence_camino,nsequence_const,totaltime,load_magnetization);

% Display function evaluation time
toc(starttime);
