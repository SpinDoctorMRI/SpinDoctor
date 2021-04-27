function setup = prepare_experiments(setup)
%PREPARE_EXPERIMENTS Prepare experiments.
%   The parameters are added to the input structure.
%
%   setup: struct
%
%   setup: struct


% We here assume we are working with water protons
gamma = 2.67513 * 1e-04;

% Check consistency of gradient sequences
nsequence = length(setup.gradient.sequences);
assert(isa(setup.gradient.sequences, "cell"))
for iseq = 1:nsequence
    seq = setup.gradient.sequences{iseq};
    assert(isa(seq, "Sequence"))
    if isa(seq, "CosOGSE") || isa(seq, "SinOGSE")
        assert(seq.nperiod > 0)
    end
end

% Assign b-values and q-values
namplitude = length(setup.gradient.values);
setup.gradient.bvalues = zeros(namplitude, nsequence);
setup.gradient.qvalues = zeros(namplitude, nsequence);
for iseq = 1:nsequence
    bnq = setup.gradient.sequences{iseq}.bvalue_no_q;
    switch setup.gradient.values_type
        case "g"
            setup.gradient.qvalues(:, iseq) = setup.gradient.values / gamma;
            setup.gradient.bvalues(:, iseq) = bnq * setup.gradient.qvalues(:, iseq).^2;
        case "q"
            setup.gradient.qvalues(:, iseq) = setup.gradient.values;
            setup.gradient.bvalues(:, iseq) = bnq * setup.gradient.qvalues(:, iseq).^2;
        case "b"
            setup.gradient.bvalues(:, iseq) = setup.gradient.values;
            setup.gradient.qvalues(:, iseq) = sqrt(setup.gradient.bvalues(:, iseq) / bnq);
        otherwise
            error("The values must be of type ""g"", ""q"" or ""b"".");
    end
end

% Create gradient directions
ndirection = setup.gradient.ndirection;
if ndirection == 1
    % Create structure for storing the one diffusion-encoding direction
    setup.gradient.directions = create_directions_onedir(setup.gradient.direction);
else
    % Obtain multiple diffusion-encoding directions uniformly distributed
    % in the unit circle or the unit sphere
    setup.gradient.directions = create_directions(ndirection, ...
        setup.gradient.flat_dirs, setup.gradient.remove_opposite);
end
assert(isequal(setup.gradient.directions.indices, 1:length(setup.gradient.directions.indices)), ...
    "directions.indices must be montonically increasing integers, starting from 1");

% Check BTPDE experiment
if isfield(setup, "btpde")
    if ~isfield(setup.btpde, "ode_solver")
        % Set default
        setup.btpde.ode_solver = @ode15s;
    elseif ischar(setup.btpde.ode_solver) || isstring(setup.btpde.ode_solver)
        setup.btpde.ode_solver = str2func(setup.btpde.ode_solver);
    elseif ~isa(setup.btpde.ode_solver, "function_handle")
        error("The BTPDE ODE solver must be a function handle or a string");
    end
end

% Check HADC experiment
if isfield(setup, "hadc")
    if ~isfield(setup.hadc, "ode_solver")
        % Set default
        setup.hadc.ode_solver = @ode15s;
    elseif ischar(setup.hadc.ode_solver) || isstring(setup.hadc.ode_solver)
        setup.hadc.ode_solver = str2func(setup.hadc.ode_solver);
    elseif ~isa(setup.hadc.ode_solver, "function_handle")
        error("The HADC ODE solver must be a function handle or a string");
    end
end

% Check MF experiment
if isfield(setup, "mf")
    assert(setup.mf.length_scale >= 0)
    assert(setup.mf.neig_max >= 1)
    assert(setup.mf.ninterval >= 1)
end

% Check analytical experiment
if isfield(setup, "analytical")
    assert(setup.analytical.length_scale >= 0)
    assert(setup.analytical.eigstep > 0)
end

% Check Karger experiment
if isfield(setup, "karger")
    if ~isfield(setup.karger, "ode_solver")
        % Set default
        setup.karger.ode_solver = @ode43;
    elseif ischar(setup.karger.ode_solver) || isstring(setup.karger.ode_solver)
        setup.karger.ode_solver = str2func(setup.karger.ode_solver);
    elseif ~isa(setup.karger.ode_solver, "function_handle")
        error("The Karger ODE solver must be a function handle or a string");
    end
end
