function setup = prepare_experiments(setup)
%PREPARE_EXPERIMENTS Prepare experiments.
%   The parameters are added to the input structure.
%
%   setup: struct
%
%   setup: struct


% Check consistency of gradient sequences
assert(isa(setup.gradient.sequences, "cell"), ...
        "Gradient sequences must use cell type data container.")
setup.gradient.sequences = setup.gradient.sequences(:)';
nsequence = length(setup.gradient.sequences);
for iseq = 1:nsequence
    seq = setup.gradient.sequences{iseq};
    assert(isa(seq, "Sequence"), ...
        "Gradient sequence must be an instance of Sequence class.")
end

% Assign b-values, q-values and g-values
namplitude = length(setup.gradient.values);
setup.gradient.bvalues = zeros(namplitude, nsequence);
setup.gradient.qvalues = zeros(namplitude, nsequence);
setup.gradient.gvalues = zeros(namplitude, nsequence);
for iseq = 1:nsequence
    bnq = setup.gradient.sequences{iseq}.bvalue_no_q;
    switch setup.gradient.values_type
        case "g"
            % commonly used magnetic field gradient unit is mT/m
            setup.gradient.qvalues(:, iseq) = setup.gradient.values * setup.gamma / 1e6;
            setup.gradient.bvalues(:, iseq) = bnq * setup.gradient.qvalues(:, iseq).^2;
        case "q"
            % q = gamma*g, different from the commonly used definition (q=gamma*g*delta)
            setup.gradient.qvalues(:, iseq) = setup.gradient.values;
            setup.gradient.bvalues(:, iseq) = bnq * setup.gradient.qvalues(:, iseq).^2;
        case "b"
            setup.gradient.bvalues(:, iseq) = setup.gradient.values;
            setup.gradient.qvalues(:, iseq) = sqrt(setup.gradient.bvalues(:, iseq) / bnq);
        otherwise
            error("The values must be of type ""g"", ""q"" or ""b"".");
    end
end
setup.gradient.gvalues = setup.gradient.qvalues * 1e6 / setup.gamma;

% Normalize gradient directions
setup.gradient.directions ...
    = setup.gradient.directions ./ vecnorm(setup.gradient.directions);

% Check BTPDE experiment
if isfield(setup, "btpde")
    setup.btpde = check_btpde(setup.btpde);
end

% Check BTPDE_MIDPOINT experiment
if isfield(setup, "btpde_midpoint")
    setup.btpde_midpoint = check_btpde_midpoint(setup.btpde_midpoint);
end

% Check HADC experiment
if isfield(setup, "hadc")
    setup.hadc = check_hadc(setup.hadc);
end

% Check MF experiment
if isfield(setup, "mf")
    setup.mf = check_mf(setup.mf);
end

% Check analytical experiment
if isfield(setup, "analytical")
    assert(setup.analytical.length_scale >= 0, ...
        "analytical: length scale is negative.")
    assert(setup.analytical.eigstep > 0, ...
        "analytical: eigstep is negative or zero.")
end

% Check Karger experiment
if isfield(setup, "karger")
    if ~isfield(setup.karger, "ode_solver")
        % Set default
        setup.karger.ode_solver = @ode45;
    elseif ischar(setup.karger.ode_solver) || isstring(setup.karger.ode_solver)
        setup.karger.ode_solver = str2func(setup.karger.ode_solver);
    elseif ~isa(setup.karger.ode_solver, "function_handle")
        error("karger: ODE solver must be a function handle or a string");
    end
end
end

function btpde = check_btpde(btpde)
    if ~isfield(btpde, "ode_solver")
        % Set default
        btpde.ode_solver = @ode15s;
    elseif ischar(btpde.ode_solver) || isstring(btpde.ode_solver)
        btpde.ode_solver = str2func(btpde.ode_solver);
    elseif ~isa(btpde.ode_solver, "function_handle")
        error("btpde: ODE solver must be a function handle or a string");
    end

    if isfield(btpde, "reltol")
        assert(btpde.reltol > 0, "btpde: reltol is negative or zero.");
    else
        btpde.reltol = 1e-4;
    end

    if isfield(btpde, "abstol")
        assert(btpde.abstol > 0, "btpde: abstol is negative or zero.");
    else
        btpde.abstol = 1e-6;
    end

    if ~isfield(btpde, "rerun")
        btpde.rerun = false;
    end
end

function btpde_mp = check_btpde_midpoint(btpde_mp)
    assert(btpde_mp.implicitness > 0, ...
        "btpde_midpoint: implicitness is negative or zero.");
    assert(btpde_mp.timestep > 0, ...
        "btpde_midpoint: timestep is negative or zero.");
    if ~isfield(btpde_mp, "rerun")
        btpde_mp.rerun = false;
    end
end

function hadc = check_hadc(hadc)
    if ~isfield(hadc, "ode_solver")
        % Set default
        hadc.ode_solver = @ode15s;
    elseif ischar(hadc.ode_solver) || isstring(hadc.ode_solver)
        hadc.ode_solver = str2func(hadc.ode_solver);
    elseif ~isa(hadc.ode_solver, "function_handle")
        error("hadc: ODE solver must be a function handle or a string");
    end

    if isfield(hadc, "reltol")
        assert(hadc.reltol > 0, "hadc: reltol is negative or zero.");
    else
        hadc.reltol = 1e-4;
    end

    if isfield(hadc, "abstol")
        assert(hadc.abstol > 0, "hadc: abstol is negative or zero.");
    else
        hadc.abstol = 1e-6;
    end

    if ~isfield(hadc, "rerun")
        hadc.rerun = false;
    end
end

function mf = check_mf(mf)
    mf.neig_max = round(mf.neig_max);
    assert(mf.neig_max >= 1, ...
        "matrix formalism: maximum number of eigenvalues must be greater than 1.");

    if isfield(mf, 'ninterval')
        mf.ninterval = round(mf.ninterval);
        assert(mf.ninterval >= 1, ...
            "matrix formalism: number of intervals must be greater than 1.");
    else
        mf.ninterval = 500;
    end

    if isfield(mf, 'length_scale')
        assert(mf.length_scale >= 0, ...
            "matrix formalism: length scale is negative.");
    else
        mf.length_scale = 0;
    end

    if isfield(mf, 'maxiter')
        mf.maxiter = round(mf.maxiter);
        assert(mf.maxiter > 0, ...
            "matrix formalism: maximum number of iterations is negative.");
    end

    if isfield(mf, 'maxiterations')
        mf.maxiterations = round(mf.maxiterations);
        assert(mf.maxiterations > 0, ...
            "matrix formalism: maximum number of iterations is negative.");
    end

    if isfield(mf, 'ssdim')
        mf.ssdim = round(mf.ssdim);
        assert(mf.ssdim > 0, ...
            "matrix formalism: maximum size of Krylov subspace is negative.");
    end

    if isfield(mf, 'subspacedimension')
        mf.subspacedimension = round(mf.subspacedimension);
        assert(mf.subspacedimension > 0, ...
            "matrix formalism: maximum size of Krylov subspace is negative.");
    end

    if isfield(mf, 'tolerance')
        assert(mf.tolerance > 0, ...
            "matrix formalism: convergence tolerance is negative.");
    end

    if ~isfield(mf, 'rerun');           mf.rerun = false;           end
    if ~isfield(mf, 'rerun_eigen');     mf.rerun_eigen = false;     end
end
