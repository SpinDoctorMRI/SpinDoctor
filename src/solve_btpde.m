function results = solve_btpde(femesh, setup)
%SOLVE_BTPDE Solve the Bloch-Torrey partial differential equation.
%
%   femesh: struct
%   setup: struct
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

% Extract HARDI directions
dir_points = setup.gradient.directions.points;
dir_inds = setup.gradient.directions.indices;
opposite = setup.gradient.directions.opposite;

% Extract domain parameters
diffusivity = setup.pde.diffusivity;
relaxation = setup.pde.relaxation;
initial_density = setup.pde.initial_density;
permeability = setup.pde.permeability;
boundary_markers = setup.pde.boundary_markers;

% Extract experiment parameters
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
reltol = setup.btpde.reltol;
abstol = setup.btpde.abstol;
solve_ode = setup.btpde.ode_solver;
solver_str = func2str(solve_ode);

% Sizes
ncompartment = length(setup.pde.compartments);
namplitude = size(qvalues, 1);
nsequence = length(sequences);
ndirection = setup.gradient.ndirection;
ndirection_unique = length(dir_inds);

% Number of points in each compartment
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

% Assemble finite element matrices
disp("Setting up FEM matrices");
M_cmpts = cell(1, ncompartment);
K_cmpts = cell(1, ncompartment);
Q_cmpts = cell(1, ncompartment);
Jx_cmpts = cell(1, 3);
for idim = 1:3
    Jx_cmpts{idim} = cell(1, ncompartment);
end
rho_cmpts = cell(1, ncompartment);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    facets = femesh.facets(icmpt, :);
    elements = femesh.elements{icmpt};
    [~, volumes] = get_volume_mesh(points, elements);

    % Assemble flux, stiffness and mass matrices in compartment
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
    K_cmpts{icmpt} = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));
    Q_cmpts{icmpt} = assemble_flux_matrix(points, facets, permeability);

    % Assemble moment matrices (coordinate weighted mass matrices)
    for idim = 1:3
        Jx_cmpts{idim}{icmpt} = mass_matrixP1_3D(elements', volumes', points(idim, :)');
    end

    % Add T2-relaxation term if it is finite and nonnegative. This term is added
    % directly to the stiffness matrix, as it has the effect of adding a
    % decay operator to the diffusion operator
    if isfinite(relaxation(icmpt)) && relaxation(icmpt) > 0
        fprintf("Adding T2-relaxation matrix in compartment %d of %d: %g\n", icmpt, ncompartment, relaxation(icmpt));
        K_cmpts{icmpt} = K_cmpts{icmpt} + 1 / relaxation(icmpt) * M_cmpts{icmpt};
    end

    % Create initial conditions (enforce complex values)
    rho_cmpts{icmpt} = complex(initial_density(icmpt)) * ones(npoint_cmpts(icmpt), 1);
end

% Create global mass, stiffness, flux and moment matrices (sparse)
disp("Coupling FEM matrices");
M = blkdiag(M_cmpts{:});
K = blkdiag(K_cmpts{:});
Q = couple_flux_matrix(femesh, boundary_markers, Q_cmpts);
Jx = cellfun(@(J) blkdiag(J{:}), Jx_cmpts, "UniformOutput", false);

% Global initial conditions
rho = vertcat(rho_cmpts{:});

% Set parameters for ODE solver
options_template = odeset("Mass", M, "AbsTol", abstol, "RelTol", reltol, ...
    "Vectorized", "on", "Stats", "off", "MassSingular", "no");

% Cartesian indices (for parallel looping with linear indices)
allinds = [namplitude nsequence ndirection];

% Iterate over gradient amplitudes, sequences and directions. If the Matlab
% PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
% parallel, otherwise it should work like a normal loop. If that is not the
% case, replace the `parfor` keyword by the normal `for` keyword.
parfor iall = 1:prod(allinds)

    % Measure iteration time
    itertime = tic;

    % Extract Cartesian indices
    [iamp, iseq, idir] = ind2sub(allinds, iall);

    % Do not bother solving if opposite direction has already been computed
    if idir <= ndirection_unique

        % Extract parameters for iteration
        q = qvalues(iamp, iseq);
        b = bvalues(iamp, iseq);
        seq = sequences{iseq};
        g = dir_points(:, idir);

        % Get intervals based on the properties of the time profile
        [timelist, interval_str, timeprofile_str] = seq.intervals;

        % Number of intervals
        ninterval = length(timelist) - 1;

        % Assemble gradient direction dependent finite element matrix
        J = g(1) * Jx{1} + g(2) * Jx{2} + g(3) * Jx{3};

        % Initial magnetization
        mag = rho;

        % Base information about current iteration
        iteration_str = sprintf("Solving BTPDE of size %d using %s\n" ...
            + "  Direction %d of %d: g = [%.2f; %.2f; %.2f]\n" ...
            + "  Sequence  %d of %d: f = %s\n" ...
            + "  Amplitude %d of %d: q = %g, b = %g", ...
            sum(npoint_cmpts), solver_str, ...
            idir, ndirection, g, ...
            iseq, nsequence, seq2str(seq), ...
            iamp, namplitude, q, b);

        % Solve for each interval consecutively
        for iint = 1:ninterval

            % Add a third point to the interval, so that the ODE solver does not
            % store the magnetization for all time steps during the solve. If
            % there were only two points in the interval, the ODE solver would
            % store all time steps. This would require a lot of memory,
            % especially during parfor iterations
            interval_midpoint = (timelist(iint) + timelist(iint + 1)) / 2;
            time_list_interval = [timelist(iint), interval_midpoint, timelist(iint + 1)];

            % Display state of iterations
            fprintf("%s\n" ...
                + "  Interval %d of %d: I = %s, %s\n", ...
                iteration_str, ...
                iint, ninterval, interval_str(iint), timeprofile_str(iint));

            % Create new ODE functions on given interval
            [ode_function, Jacobian] = btpde_functions_interval(K, Q, J, q, seq, interval_midpoint);

            % Update options with new Jacobian, which is either a
            % function handle or a constant matrix, depending on the
            % time profile
            options = odeset(options_template, "Jacobian", Jacobian);

            % Solve ODE on domain, starting from the magnetization at
            % the end of the previous interval (mag)
            [~, y] = solve_ode(ode_function, time_list_interval, mag, options);

            % Magnetization at end of interval
            mag = y(end, :).';
        end

        % Split global solution into compartments
        mag = mat2cell(mag, npoint_cmpts).';
        % magnetization(:, iall) = mag;  % does not work in parfor before R2020b
        for icmpt = 1:ncompartment
            magnetization{icmpt, iall} = mag{icmpt};
        end
        % [magnetization{:, iall}] = deal(mag{:});
        signal(:, iall) = cellfun(@(M, y) sum(M * y, 1), M_cmpts, mag);

    end % unique direction case

    % Store timing
    itertimes(iall) = toc(itertime);

end % iterations

% Copy results to opposite directions
for idir = 1:ndirection_unique
    if ~isempty(opposite{idir})
        fprintf("Copying complex conjugate results from direction %d to opposite direction (%d)\n", idir, opposite{idir});
        magnetization(:, :, :, opposite{idir}) = cellfun(@conj, magnetization(:, :, :, idir), ...
            "UniformOutput", false);
        signal(:, :, :, opposite{idir}) = conj(signal(:, :, :, idir));
    end
end

% Total magnetization (sum over compartments)
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.magnetization = magnetization;
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime);

% Display function evaluation time
toc(starttime);
