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
    assert(isa(seq, "AbsSequence"), ...
        "Gradient sequence must be an instance of Sequence class.")
end
const_seq_ind = cellfun(@(s) ~isa(s,'SequenceCamino'),setup.gradient.sequences,'UniformOutput',true);
nsequence_const = sum(const_seq_ind);

% Assign b-values, q-values and g-values
if nsequence_const > 0 
    const_seq= setup.gradient.sequences(const_seq_ind);
    namplitude = length(setup.gradient.values);
    setup.gradient.bvalues = zeros(namplitude, nsequence_const);
    setup.gradient.qvalues = zeros(namplitude, nsequence_const);
    setup.gradient.gvalues = zeros(namplitude, nsequence_const);
    for iseq = 1:nsequence_const
        bnq = const_seq{iseq}.bvalue_no_q;
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
else
    setup.gradient.qvalues =[];
    setup.gradient.bvalues =[];
    setup.gradient.gvalues =[];
    setup.gradient.directions =[];
end
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
    if numel(unique(setup.pde.initial_density)) ~= 1
        warning("mf: matrix formalism method is inaccurate " + ...
            "for nonequilibrium spin density distribution.");
    end
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

    if isfield(mf, 'single')
        if mf.single
            mf.single = true;
            msg = "matrix formalism: using single precision.";
            disp(msg)
        else
            mf.single = false;
        end
    else
        mf.single = false;
    end

    if isfield(mf, 'gpu')
        if ~mf.gpu
            mf.gpu = false;
        else
            if ~license('test', 'Distrib_Computing_Toolbox')
                mf.gpu = false;
                msg = "matrix formalism: Parallel Computing Toolbox is not available, use CPU.";
                warning(msg)
            else
                availableGPUs = gpuDeviceCount("available");

                if availableGPUs==0
                    mf.gpu = false;
                    warning("matrix formalism: no GPU is available, use CPU.")
                else
                    if isnumeric(mf.gpu) && (numel(mf.gpu) > 1)
                        mf.gpu = unique(reshape(mf.gpu, 1, []));
                        if any(mf.gpu<=0, 'all')
                            msg = "matrix formalism: GPU index starts from 1.";
                            error(msg)
                        end
                        if any(mf.gpu>availableGPUs, 'all')
                            msg = "matrix formalism: GPU index exceeds the number of GPU devices.";
                            error(msg)
                        end
                    elseif isnumeric(mf.gpu) && (numel(mf.gpu) == 1)
                        mf.gpu = round(mf.gpu);
                        if mf.gpu <= 0
                            msg = "matrix formalism: GPU index starts from 1.";
                            error(msg)
                        end
                        if mf.gpu > availableGPUs
                            msg = "matrix formalism: GPU index exceeds the number of GPU devices.";
                            error(msg)
                        end
                    else
                        mf.gpu = true;
                    end
                end
            end
        end
    else
        mf.gpu = false;
    end % isfield

    if isfield(mf, 'length_scale')
        assert(mf.length_scale >= 0, ...
            "matrix formalism: length scale is negative.");
    else
        mf.length_scale = 0;
    end

    if isfield(mf, 'ninterval')
        mf.ninterval = round(mf.ninterval);
        assert(mf.ninterval >= 1, ...
            "matrix formalism: number of intervals must be greater than 1.");
    else
        mf.ninterval = 500;
    end

    if ~isfield(mf, 'rerun');           mf.rerun = false;           end
    if ~isfield(mf, 'rerun_eigen');     mf.rerun_eigen = false;     end
    if ~isfield(mf, 'surf_relaxation'); mf.surf_relaxation = false; end
    
    if isinf(mf.neig_max)
        % Infinite neig_max triggers eig instead of eigs
        mf = rmfields(mf, {'eigs'});
    elseif isfield(mf, 'eigs')
        % Check Matlab's eigs settings
        % more info https://www.mathworks.com/help/matlab/ref/eigs.html
        if isfield(mf.eigs, 'sigma')
            if isstring(mf.eigs.sigma) || ischar(mf.eigs.sigma)
                mf.eigs.sigma = string(mf.eigs.sigma);
            elseif ~isnumeric(mf.eigs.sigma)
                error("mf.eigs.sigma must be a string or a scalar.");
            end
        else
            mf.eigs.sigma = -1e-8;
        end

        % alias of maxiterations
        if isfield(mf.eigs, 'maxiter')
            mf.eigs.maxiter = round(mf.eigs.maxiter);
            assert(mf.eigs.maxiter > 0, ...
                "matrix formalism: maximum number of iterations is negative.");
        end

        if isfield(mf.eigs, 'maxiterations')
            mf.eigs.maxiterations = round(mf.eigs.maxiterations);
            assert(mf.eigs.maxiterations > 0, ...
                "matrix formalism: maximum number of iterations is negative.");
        end

        % alias of subspacedimensions
        if isfield(mf.eigs, 'ssdim')
            mf.eigs.ssdim = round(mf.eigs.ssdim);
            assert(mf.eigs.ssdim > 0, ...
                "matrix formalism: maximum size of Krylov subspace is negative.");
        end

        if isfield(mf.eigs, 'subspacedimension')
            mf.eigs.subspacedimension = round(mf.eigs.subspacedimension);
            assert(mf.eigs.subspacedimension > 0, ...
                "matrix formalism: maximum size of Krylov subspace is negative.");
        end

        if isfield(mf.eigs, 'tolerance')
            assert(mf.eigs.tolerance > 0, ...
                "matrix formalism: convergence tolerance is negative.");
        end
    else
        mf.eigs.sigma = -1e-8;
    end
end
