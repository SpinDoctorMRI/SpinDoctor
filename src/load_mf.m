function results = load_mf(setup, savepath, load_magnetization)
%LOAD_MF Compute the solution to the BTPDE using Matrix Formalism.
%
%   femesh: struct
%   setup: struct
%   lap_eig: struct with fields
%       values: double(neig, 1)
%       funcs: double(npoint, neig)
%   savepath: path string
%   save_magnetization (optinal): boolean 
%
%   results: struct with fields
%       signal: double(ncompartment, nsequence, namplitude, ndirection)
%       signal_allcmpts: double(nsequence, namplitude, ndirection)
%       ctime: double(ncompartment, nsequence, namplitude, ndirection)
%       magnetization: cell(ncompartment, nsequence, namplitude, ndirection)


% Measure time of function evaluation
starttime = tic;

% Provide default value
if nargin < nargin(@load_btpde)
    load_magnetization = true;
end

% Extract experiment parameters
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;
ninterval = setup.mf.ninterval;

% Sizes
ncompartment = setup.ncompartment;
namplitude = setup.namplitude;
nsequence = setup.nsequence;
ndirection = setup.ndirection;

% Folder for saving
savepath = sprintf("%s/ninterval%d", savepath, ninterval);

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(ncompartment, namplitude, nsequence, ndirection);

% Cartesian indices (for parallel looping with linear indices)
allinds = [namplitude nsequence ndirection];

% save magnetization to parfor_mag to avoid I/O error using parfor
parfor_mag = cell(allinds);
parfor iall = 1:prod(allinds)
    % Extract indices
    [iamp, iseq, idir] = ind2sub(allinds, iall);

    % Experiment parameters
    seq = sequences{iseq};
    b = bvalues(iamp, iseq);
    ug = directions(:, idir);

    % File name for saving or loading iteration results
    filename = sprintf("%s/%s.mat", savepath, seq.string(true));

    % Load results
    fprintf("Load mf %d/%d.\n", iall, prod(allinds));
    mfile = matfile(filename, "Writable", false);
    data = mfile.(gradient_fieldstring(ug, b));

    signal(:, iall) = data.signal;
    itertimes(:, iall) = data.itertimes;
    if load_magnetization
        parfor_mag{iall} = data.magnetization;
    end
end % parfor iterations

if load_magnetization
    for iall = 1:prod(allinds)
        for icmpt = 1:ncompartment
            magnetization{icmpt, iall} = parfor_mag{iall}{icmpt};
        end
    end
end

% Compute total signal (sum over compartments)
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
