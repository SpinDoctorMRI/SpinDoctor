function experiment = prepare_experiments(experiment)
%PREPARE_EXPERIMENTS Prepare experiments.
%   The parameters are added to the input structure.
%
%   experiment: struct
%
%   experiment: struct


% We here assume we are working with water protons
gamma = 2.67513 * 1e-04;

% Check consistency of gradient sequences
nsequence = length(experiment.sequences);
assert(isa(experiment.sequences, "cell"))
for iseq = 1:nsequence
    seq = experiment.sequences{iseq};
    assert(isa(seq, "Sequence"))
    if isa(seq, "CosOGSE") || isa(seq, "SinOGSE")
        assert(seq.nperiod > 0)
    end
end

% Assign b-values and q-values
namplitude = length(experiment.values);
experiment.bvalues = zeros(namplitude, nsequence);
experiment.qvalues = zeros(namplitude, nsequence);
for iseq = 1:nsequence
    bnq = experiment.sequences{iseq}.bvalue_no_q;
    switch experiment.values_type
        case "g"
            experiment.qvalues(:, iseq) = experiment.values / gamma;
            experiment.bvalues(:, iseq) = bnq * experiment.qvalues(:, iseq).^2;
        case "q"
            experiment.qvalues(:, iseq) = experiment.values;
            experiment.bvalues(:, iseq) = bnq * experiment.qvalues(:, iseq).^2;
        case "b"
            experiment.bvalues(:, iseq) = experiment.values;
            experiment.qvalues(:, iseq) = sqrt(experiment.bvalues(:, iseq) / bnq);
        otherwise
            error("The experiment values must be of type ""g"", ""q"" or ""b"".");
    end
end

% Create gradient directions
ndirection = experiment.ndirection;
if ndirection == 1
    % Create structure for storing the one diffusion-encoding direction
    experiment.directions = create_directions_onedir(experiment.direction);
else
    % Obtain multiple diffusion-encoding directions uniformly distributed
    % in the unit circle or the unit sphere
    experiment.directions = create_directions(ndirection, ...
        experiment.flat_dirs, experiment.remove_opposite);
end
assert(isequal(experiment.directions.indices, 1:length(experiment.directions.indices)), ...
    "directions.indices must be montonically increasing integers, starting from 1");

% Check BTPDE experiment
if isfield(experiment, "btpde")
    if ~isfield(experiment.btpde, "ode_solver")
        % Set default
        experiment.btpde.("ode_solver") = @ode15s;
    elseif ischar(experiment.btpde.ode_solver)
        experiment.btpde.ode_solver = str2func(experiment.btpde.ode_solver);
    elseif ~isa(experiment.btpde.ode_solver, "function_handle")
        error("The BTPDE ODE solver must be a function handle or a character vector");
    end
end

% Check HADC experiment
if isfield(experiment, "hadc")
    if ~isfield(experiment.hadc, "ode_solver")
        % Set default
        experiment.btpde.("ode_solver") = @ode15s;
    elseif ischar(experiment.hadc.ode_solver)
        experiment.hadc.ode_solver = str2func(experiment.hadc.ode_solver);
    elseif ~isa(experiment.hadc.ode_solver, "function_handle")
        error("The HADC ODE solver must be a function handle or a character vector");
    end
end

% Check MF experiment
if isfield(experiment, "mf")
    assert(experiment.mf.length_scale >= 0)
    assert(experiment.mf.neig_max >= 1)
    assert(experiment.mf.ninterval >= 1)
end

% Check multilayer experiment
if isfield(experiment, "multilayer")
    assert(experiment.multilayer.length_scale >= 0)
    assert(experiment.multilayer.eigstep > 0)
end
