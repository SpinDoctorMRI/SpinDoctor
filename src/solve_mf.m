function results = solve_mf(femesh, setup, lap_eig)
%SOLVE_MF Compute the solution to the BTPDE using Matrix Formalism.
%
%   femesh: struct
%   setup: struct
%   lap_eig: struct with fields
%       values: cell(1, ndomain)
%       funcs: cell(1, ndomain)
%       moments: cell(1, ndomain)
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
dir_points = setup.gradient.directions.points;
dir_inds = setup.gradient.directions.indices;
opposite = setup.gradient.directions.opposite;
eigfuncs = lap_eig.funcs;
moments = lap_eig.moments;
T2 = lap_eig.massrelax;
L = diag(lap_eig.values);

% Sizes
ncompartment = femesh.ncompartment;
namplitude = length(setup.gradient.values);
nsequence = length(sequences);
ndirection = setup.gradient.ndirection;
ndirection_unique = length(dir_inds);
neig = length(lap_eig.values);
ninterval = setup.mf.ninterval;

% Number of points in each compartment
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

% Assemble mass matrix in each compartment (for spatial integration)
disp("Assembling mass matrices");
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

% Global matrix and density
M = blkdiag(M_cmpts{:});
rho = vertcat(rho_cmpts{:});

% Coefficients of initial spin density in Laplace eigenfunction basis
nu0 = lap_eig.funcs' * (M * rho);

% Cartesian indices (for parallel looping with linear indices)
allinds = [namplitude nsequence ndirection];

% Iterate over gradient amplitudes, sequences and directions. If the Matlab
% PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
% parallel, otherwise it should work like a normal loop. If that is not the
% case, replace the `parfor` keyword by the normal `for` keyword.
parfor iall = 1:prod(allinds)

    % Measure iteration time
    itertime = tic;

    % Extract indices
    [iamp, iseq, idir] = ind2sub(allinds, iall);

    % Do not bother solving if opposite direction has already been computed
    if idir <= ndirection_unique

        % Experiment parameters
        q = qvalues(iamp, iseq);
        b = bvalues(iamp, iseq);
        seq = sequences{iseq};
        g = dir_points(:, idir);

        % Create time intervals for time profile approximation
        time = linspace(0, seq.echotime, ninterval + 1);

        % Display state of iterations
        fprintf("Computing MF magnetization using %d eigenvalues\n" ...
            + "  Direction %d of %d: g = [%.2f; %.2f; %.2f]\n" ...
            + "  Sequence  %d of %d: f = %s\n" ...
            + "  Amplitude %d of %d: q = %g, b = %g\n", ...
            neig, ...
            idir, ndirection, g, ...
            iseq, nsequence, seq2str(seq), ...
            iamp, namplitude, q, b);

        % Components of BT operator matrix
        A = sum(moments .* shiftdim(g, -2), 3);

        % Compute final magnetization (in Laplace basis)
        % If PGSE, use three intervals, otherwise many intervals
        if isa(seq, "PGSE")
            % Constant BT operator in Laplace basis
            K = L + T2 + 1i * q * A;
            edK = expm(-seq.delta * K);
            edL = expm(-(seq.Delta - seq.delta) * (L + T2));

            % Laplace coefficients of final magnetization
            nu = edK' * (edL * (edK * nu0));
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
                nu = expm(-dt * K(ft)) * nu;
            end
        end

        % Final magnetization coefficients in finite element nodal basis
        mag = eigfuncs * nu;

        % Split magnetization into compartments
        mag = mat2cell(mag, npoint_cmpts);
        % magnetization(:, iall) = mag;  % does not work in parfor before R2020b
        for icmpt = 1:ncompartment
            magnetization{icmpt, iall} = mag{icmpt};
        end
        signal(:, iall) = cellfun(@(M, m) sum(M * m), M_cmpts', mag);

    end % unique direction case

    % Store timing
    itertimes(iall) = toc(itertime);

end % iterations


% Copy result to opposite directions
for idir = dir_inds
    if ~isempty(opposite{idir})
        fprintf("Copying complex conjugate results from direction %d to opposite direction (%d)\n", idir, opposite{idir});
        magnetization(:, :, :, opposite{idir}) = cellfun(@conj, magnetization(:, :, :, idir), ...
            "UniformOutput", false);
        signal(:, :, :, opposite{idir}) = conj(signal(:, :, :, idir));
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
