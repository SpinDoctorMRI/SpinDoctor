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


% Measure function evaluation time
starttime = tic;

% Provide default value
if nargin < nargin(@load_btpde)
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
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

inds = [namplitude ndirection];
% Iterate over gradient amplitudes, sequences and directions.
for iseq = 1:nsequence
    seq = sequences{iseq};
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

        signal(:, iamp, iseq, idir) = data.signal;
        itertimes(iamp, iseq, idir) = data.itertimes;
        if load_magnetization
            magnetization(:, iamp, iseq, idir) = data.magnetization;
        end
    end
end

% Total magnetization (sum over compartments)
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = sum(itertimes, "all");
if load_magnetization
    results.magnetization = magnetization;
    results.magnetization_avg = average_magnetization(magnetization);
end

% Display function evaluation time
toc(starttime);
