function results = load_btpde_midpoint(setup, savepath, load_magnetization)
%LOAD_BTPDE_MIDPOINT Load the results saved by SOLVE_BTPDE_MIDPOINT.
%
%   LOAD_BTDPE(SETUP, SAVEPATH) loads the results of each iteration from
%   "<SAVEPATH>/BTPDE_MIDPOINT_<SOLVEROPTIONS>/<ITERATIONINFO>.MAT".
%
%   LOAD_BTDPE(SETUP, SAVEPATH, LOAD_MAGNETIZATION) also omits loading
%   the magnetization field if LOAD_MAGNETIZATION is set to FALSE.
%
%   setup: struct
%   savepath (optional): string
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
values = setup.gradient.values;
amptype = setup.gradient.values_type;
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;
theta = setup.btpde_midpoint.implicitness;
dt = setup.btpde_midpoint.timestep;

% Sizes
ncompartment = length(setup.pde.initial_density);
namplitude = size(qvalues, 1);
nsequence = length(sequences);
ndirection = size(directions, 2);

% Folder for saving
savepath = sprintf( ...
    "%s/btpde_midpoint_theta%g_dt%g_magnetization%d", ...
    savepath, theta, dt, load_magnetization ...
);

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

% Cartesian indices (for parallel looping with linear indices)
allinds = [namplitude nsequence ndirection];

% Iterate over gradient amplitudes, sequences and directions. If the Matlab
% PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
% parallel, otherwise it should work like a normal loop. If that is not the
% case, replace the `parfor` keyword by the normal `for` keyword.
parfor iall = 1:prod(allinds)

    % Extract Cartesian indices
    [iamp, iseq, idir] = ind2sub(allinds, iall);

    % Extract iteration inputs
    amp = values(iamp);
    q = qvalues(iamp, iseq);
    b = bvalues(iamp, iseq);
    seq = sequences{iseq};
    g = directions(:, idir);
    
    % File name for saving or loading iteration results
    filename = sprintf("%s/%s.mat", savepath, gradient_string(amp, amptype, seq, g));
    
    % Load results
    fprintf("Load %s\n", filename);
    mfile = matfile(filename, "Writable", false);
    signal(:, iall) = mfile.signal;
    itertimes(iall) = mfile.itertime;
    if load_magnetization
        for icmpt = 1:ncompartment
            magnetization{icmpt, iall} = mfile.magnetization;
        end
    end
end % iterations

% Total magnetization (sum over compartments)
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
if load_magnetization
    results.magnetization = magnetization;
end
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = sum(itertimes, "all");

% Display function evaluation time
toc(starttime);
