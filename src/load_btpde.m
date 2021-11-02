function results = load_btpde(setup, savepath, load_magnetization)
%LOAD_BTPDE Load the results saved by SOLVE_BTPDE.
%
%   LOAD_BTDPE(SETUP, SAVEPATH) loads the results of each iteration from
%   "<SAVEPATH>/<SOLVER_OPTIONS>/<ITERATION_INFO>.MAT".
%
%   LOAD_BTDPE(SETUP, SAVEPATH, LOAD_MAGNETIZATION) also omits loading
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
reltol = setup.btpde.reltol;
abstol = setup.btpde.abstol;
solve_ode = setup.btpde.ode_solver;
solver_str = func2str(solve_ode);

% Sizes
ncompartment = setup.ncompartment;
namplitude = setup.namplitude;
nsequence = setup.nsequence;
ndirection = setup.ndirection;

% Folder for saving
savepath = sprintf( ...
    "%s/%s_abstol%g_reltol%g", ...
    savepath, solver_str, abstol, reltol ...
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

% Save magnetization in temp_mag to avoid parfor IO error
temp_mag = cell(namplitude, nsequence, ndirection);
parfor iall = 1:prod(allinds)
    % Extract Cartesian indices
    [iamp, iseq, idir] = ind2sub(allinds, iall);

    % Extract iteration inputs
    seq = sequences{iseq};
    b = bvalues(iamp, iseq);
    ug = directions(:, idir);

    % File name for saving or loading iteration results
    filename = sprintf("%s/%s.mat", savepath, seq.string(true));

    % Load results
    fprintf("Load btpde %d/%d.\n", iall, prod(allinds));
    mfile = matfile(filename, "Writable", false);
    data = mfile.(gradient_fieldstring(ug, b));
    
    signal(:, iall) = data.signal;
    itertimes(iall) = data.itertimes;
    if load_magnetization
        temp_mag{iall} = data.magnetization;
    end
end % iterations

for iall = 1:prod(allinds)
    magnetization(:, iall) = temp_mag{iall};
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
