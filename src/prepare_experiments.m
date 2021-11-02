function setup = prepare_experiments(setup)
%PREPARE_EXPERIMENTS Prepare experiments.
%   The parameters are added to the input structure.
%
%   setup: struct
%
%   setup: struct


% Check consistency of gradient sequences
assert(isa(setup.gradient.sequences, "cell"))
setup.gradient.sequences = setup.gradient.sequences(:)';
nsequence = length(setup.gradient.sequences);
for iseq = 1:nsequence
    seq = setup.gradient.sequences{iseq};
    assert(isa(seq, "Sequence"))
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
    assert(setup.analytical.length_scale >= 0)
    assert(setup.analytical.eigstep > 0)
end

% Check Karger experiment
if isfield(setup, "karger")
    if ~isfield(setup.karger, "ode_solver")
        % Set default
        setup.karger.ode_solver = @ode45;
    elseif ischar(setup.karger.ode_solver) || isstring(setup.karger.ode_solver)
        setup.karger.ode_solver = str2func(setup.karger.ode_solver);
    elseif ~isa(setup.karger.ode_solver, "function_handle")
        error("The Karger ODE solver must be a function handle or a string");
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
        error("The BTPDE ODE solver must be a function handle or a string");
    end

    if isfield(btpde, "reltol")
        assert(btpde.reltol > 0);
    else
        btpde.reltol = 1e-4;
    end

    if isfield(btpde, "abstol")
        assert(btpde.abstol > 0);
    else
        btpde.abstol = 1e-6;
    end

    if ~isfield(btpde, "rerun")
        btpde.rerun = false;
    end
end

function btpde_mp = check_btpde_midpoint(btpde_mp)
    assert(btpde_mp.implicitness > 0);
    assert(btpde_mp.timestep > 0);
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
        error("The HADC ODE solver must be a function handle or a string");
    end

    if isfield(hadc, "reltol")
        assert(hadc.reltol > 0);
    else
        hadc.reltol = 1e-4;
    end

    if isfield(hadc, "abstol")
        assert(hadc.abstol > 0);
    else
        hadc.abstol = 1e-6;
    end

    if ~isfield(hadc, "rerun")
        hadc.rerun = false;
    end
end

function mf = check_mf(mf)
    assert(mf.neig_max >= 1);

    if isfield(mf, 'ninterval')
        assert(mf.ninterval >= 1);
    else
        mf.ninterval = 500;
    end

    if isfield(mf, 'length_scale')
        assert(mf.length_scale >= 0);
    else
        mf.length_scale = 0;
    end

    if ~isfield(mf, 'rerun');               mf.rerun = false;                   end
    if ~isfield(mf, 'rerun_eigen');         mf.rerun_eigen = false;             end
    if isfield(mf, 'tolerance');            assert(mf.tolerance > 0);           end
    if isfield(mf, 'maxiter');              assert(mf.maxiter > 0);             end
    if isfield(mf, 'maxiterations');        assert(mf.maxiterations > 0);       end
    if isfield(mf, 'ssdim');                assert(mf.ssdim > 0);               end
    if isfield(mf, 'subspacedimension');    assert(mf.subspacedimension > 0);   end
end
