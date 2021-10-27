function results = solve_mf(femesh, setup, lap_eig)
%SOLVE_MF Compute the solution to the BTPDE using Matrix Formalism.
%
%   femesh: struct
%   setup: struct
%   lap_eig: struct with fields
%       values: double(neig, 1)
%       funcs: double(npoint, neig)
%       moments: double(neig, neig, 3)
%       massrelax: double(neig, neig)
%
%   results: struct with fields
%       signal: double(ncompartment, nsequence, namplitude, ndirection)
%       signal_allcmpts: double(nsequence, namplitude, ndirection)
%       ctime: double(nsequence, namplitude, ndirection)
%       magnetization: cell(nsequence, namplitude, ncompartment, ndirection)


% Measure time of function evaluation
starttime = tic;

% Extract parameters
initial_density = setup.pde.initial_density;
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;

% Sizes
ncompartment = femesh.ncompartment;
namplitude = length(setup.gradient.values);
nsequence = length(sequences);
ndirection = size(directions, 2);
ninterval = setup.mf.ninterval;

% Cartesian indices (for parallel looping with linear indices)
allinds = [namplitude nsequence ndirection];

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

% Assemble mass matrix in each compartment (for spatial integration)
M_cmpts = cell(1, ncompartment);
rho_cmpts = cell(1, ncompartment);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [~, volumes] = get_volume_mesh(points, elements);

    % Mass matrix
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');

    % Create initial conditions (enforce complex values)
    rho_cmpts{icmpt} = complex(initial_density(icmpt)) * ones(size(femesh.points{icmpt}, 2), 1);
end

if length(lap_eig) == 1    % One compartment or some compartments are connected by permeable interfaces
    eigfuncs = lap_eig.funcs;
    moments = lap_eig.moments;
    T2 = lap_eig.massrelax;
    L = diag(lap_eig.values);
    neig = length(lap_eig.values);

    % Number of points in each compartment
    npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

    % Global matrix and density
    M = blkdiag(M_cmpts{:});
    rho = vertcat(rho_cmpts{:});
    % Coefficients of initial spin density in Laplace eigenfunction basis
    nu0 = eigfuncs' * (M * rho);

    % Iterate over gradient amplitudes, sequences and directions. If the Matlab
    % PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
    % parallel, otherwise it should work like a normal loop. If that is not the
    % case, replace the `parfor` keyword by the normal `for` keyword.
    parfor iall = 1:prod(allinds)

        % Measure iteration time
        itertime = tic;

        % Extract indices
        [iamp, iseq, idir] = ind2sub(allinds, iall);

        % Experiment parameters
        q = qvalues(iamp, iseq);
        b = bvalues(iamp, iseq);
        seq = sequences{iseq};
        g = directions(:, idir);

        % Create time intervals for time profile approximation
        time = linspace(0, seq.echotime, ninterval + 1);

        % Display state of iterations
        fprintf("Computing MF magnetization using %d eigenvalues\n" ...
            + "  Direction %d of %d: g = [%.2f; %.2f; %.2f]\n" ...
            + "  Sequence  %d of %d: f = %s\n" ...
            + "  Amplitude %d of %d: q = %g, b = %g\n", ...
            neig, ...
            idir, ndirection, g, ...
            iseq, nsequence, seq, ...
            iamp, namplitude, q, b);

        % Components of BT operator matrix
        A = sum(moments .* shiftdim(g, -2), 3);

        % Compute final magnetization (in Laplace basis)
        % If PGSE, use three intervals; if DoublePGSE, use seven intervals; otherwise many intervals.
        if isa(seq, "PGSE")
            % Constant BT operator in Laplace basis
            K = L + T2 + 1i * q * A;

            % Laplace coefficients of final magnetization
            nu = expmv(-seq.delta, K, nu0);
            nu = expm(-(seq.Delta - seq.delta) * (L + T2)) * nu;
            nu = expmv(-seq.delta, K', nu);

            % edK = expm(-seq.delta * K);
            % edL = expm(-(seq.Delta - seq.delta) * (L + T2));
            % nu = edK' * (edL * (edK * nu0));
        elseif isa(seq, "DoublePGSE")
            % Constant BT operator in Laplace basis
            K = L + T2 + 1i * q * A;

            % Laplace coefficients of final magnetization
            edK = expm(-seq.delta * K);
            edL = expm(-(seq.Delta - seq.delta) * (L + T2));
            etL = expm(-seq.tpause * (L + T2));
            nu = edK' * (edL * (edK * (etL * (edK' * (edL * (edK * nu0))))));

            % nu = expmv(-seq.delta, K, nu0);
            % nu = expm(-(seq.Delta - seq.delta) * (L + T2)) * nu;
            % nu = expmv(-seq.delta, K', nu);
            % nu = expm(-seq.tpause * (L + T2)) * nu;
            % nu = expmv(-seq.delta, K, nu);
            % nu = expm(-(seq.Delta - seq.delta) * (L + T2)) * nu;
            % nu = expmv(-seq.delta, K', nu);
        else
            % BT operator in Laplace basis for a given time profile value
            K = @(ft) L + T2 + 1i * q * ft * A;

            % Transform Laplace coefficients using piecewise constant
            % approximation of time profile
            nu = nu0;
            for i = 1:ninterval
                % Time step and time profile on given interval
                dt = time(i + 1) - time(i);
                ft = (seq.call(time(i + 1)) + seq.call(time(i))) / 2;

                % Laplace coefficients of magnetization at end of interval
                nu = expmv(-dt, K(ft), nu);
                % nu = expm(-dt * K(ft)) * nu;
            end
        end

        % Final magnetization coefficients in finite element nodal basis
        mag = eigfuncs * nu;
        
        % Split magnetization into compartments
        mag = mat2cell(mag, npoint_cmpts);
        for icmpt = 1:ncompartment
            magnetization{icmpt, iall} = mag{icmpt};
        end
        signal(:, iall) = cellfun(@(M, m) sum(M * m), M_cmpts', mag);

        % Store timing
        itertimes(iall) = toc(itertime);
    end % parfor iterations

else    % All compartments are uncorrelated
    for icmpt = 1:ncompartment
        % Extract eigenvalues, eigenfunctions, moments and T2-relaxation matrix
        eigfuncs = lap_eig(icmpt).funcs;
        moments = lap_eig(icmpt).moments;
        T2 = lap_eig(icmpt).massrelax;
        L = diag(lap_eig(icmpt).values);
        neig = length(lap_eig(icmpt).values);
        
        % Display state of iterations
        fprintf("Computing MF magnetization for compartment %d " ...
            + "using %d eigenvalues.\n", icmpt, neig);

        % Mass matrix and density
        M = M_cmpts{icmpt};
        rho = rho_cmpts{icmpt};
        % Coefficients of initial spin density in Laplace eigenfunction basis
        nu0 = eigfuncs' * (M * rho);

        parfor iall = 1:prod(allinds)

            % Measure iteration time
            itertime = tic;

            % Extract indices
            [iamp, iseq, idir] = ind2sub(allinds, iall);

            % Experiment parameters
            q = qvalues(iamp, iseq);
            b = bvalues(iamp, iseq);
            seq = sequences{iseq};
            g = directions(:, idir);

            % Create time intervals for time profile approximation
            time = linspace(0, seq.echotime, ninterval + 1);

            % Components of BT operator matrix
            A = sum(moments .* shiftdim(g, -2), 3);

            % Compute final magnetization (in Laplace basis)
            % If PGSE, use three intervals; if DoublePGSE, use seven intervals; otherwise many intervals.
            if isa(seq, "PGSE")
                % Constant BT operator in Laplace basis
                K = L + T2 + 1i * q * A;
    
                % Laplace coefficients of final magnetization
                nu = expmv(-seq.delta, K, nu0);
                nu = expm(-(seq.Delta - seq.delta) * (L + T2)) * nu;
                nu = expmv(-seq.delta, K', nu);
    
                % edK = expm(-seq.delta * K);
                % edL = expm(-(seq.Delta - seq.delta) * (L + T2));
                % nu = edK' * (edL * (edK * nu0));
            elseif isa(seq, "DoublePGSE")
                % Constant BT operator in Laplace basis
                K = L + T2 + 1i * q * A;
    
                % Laplace coefficients of final magnetization
                edK = expm(-seq.delta * K);
                edL = expm(-(seq.Delta - seq.delta) * (L + T2));
                etL = expm(-seq.tpause * (L + T2));
                nu = edK' * (edL * (edK * (etL * (edK' * (edL * (edK * nu0))))));
    
                % nu = expmv(-seq.delta, K, nu0);
                % nu = expm(-(seq.Delta - seq.delta) * (L + T2)) * nu;
                % nu = expmv(-seq.delta, K', nu);
                % nu = expm(-seq.tpause * (L + T2)) * nu;
                % nu = expmv(-seq.delta, K, nu);
                % nu = expm(-(seq.Delta - seq.delta) * (L + T2)) * nu;
                % nu = expmv(-seq.delta, K', nu);
            else
                % BT operator in Laplace basis for a given time profile value
                K = @(ft) L + T2 + 1i * q * ft * A;
    
                % Transform Laplace coefficients using piecewise constant
                % approximation of time profile
                nu = nu0;
                for i = 1:ninterval
                    % Time step and time profile on given interval
                    dt = time(i + 1) - time(i);
                    ft = (seq.call(time(i + 1)) + seq.call(time(i))) / 2;
    
                    % Laplace coefficients of magnetization at end of interval
                    nu = expmv(-dt, K(ft), nu);
                    % nu = expm(-dt * K(ft)) * nu;
                end
            end

            % Final magnetization coefficients in finite element nodal basis
            mag = eigfuncs * nu;
            magnetization{icmpt, iall} = mag;
            signal(icmpt, iall) = sum(M * mag);

            % Store timing
            itertimes(iall) = itertimes(iall) + toc(itertime);
        end % parfor iterations
    end
end

% Compute total signal (sum over compartments)
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.magnetization = magnetization;
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime);

% Display function evaluation time
toc(starttime);
