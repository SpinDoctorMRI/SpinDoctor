function results = solve_btpde_midpoint(femesh, setup)
%SOLVE_BTPDE_MIDPOINT Solve the Bloch-Torrey partial differential equation.
%   This solver uses a theta timestepping method (generalized midpoint), and
%   only works for interval-wise constant time profiles (PGSE, DoublePGSE).
%
%   The ODE is defined by
%
%       M * dy/dt = J(t) * y,
%
%   where `J(t) = J` is interval-wise constant. The time stepping is given by
%       
%       M*(y-y0)/dt = theta*J*y + (1-theta)*J*y0,
%   
%   where `y0` is the current solution and `y` is the next solution. This is
%   rewritten
%
%       (M - dt*theta*J) * y = (M + dt*(1-theta)*J) * y0.
%
%   `J` being constant on the entire interval, we can decompose
%
%       (M - dt*theta*J) = L*U
%
%   where `L` and `U` are lower and upper triangular matrices respectively.
%   Given a constant time step `dt`, this factorization only needs to be done
%   once per interval.
%   In constrast, ODE15S uses adaptive time stepping, leading to a number of
%   refactorizations. ODE15S does however manage to reuse the LU-factorization
%   quite a few steps, by not changing the step sizes to often.
%
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

% Extract domain parameters
diffusivity = setup.pde.diffusivity;
relaxation = setup.pde.relaxation;
initial_density = setup.pde.initial_density;

% Extract experiment parameters
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;
theta = setup.btpde_midpoint.implicitness;
dt = setup.btpde_midpoint.timestep;
solver_str = sprintf("theta rule (theta = %g, dt = %g)", theta, dt);

% Check that sequences are compatible
assert(all(cellfun(@(f) isa(f, "PGSE") || isa(f, "DoublePGSE"), sequences)), ...
    "Only constant sequences are currently supported.");

% Sizes
ncompartment = femesh.ncompartment;
namplitude = size(qvalues, 1);
nsequence = length(sequences);
ndirection = size(directions, 2);

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
R_cmpts = cell(1, ncompartment);
Jx_cmpts = repmat({cell(1, ncompartment)}, 1, 3);
rho_cmpts = cell(1, ncompartment);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [~, volumes] = get_volume_mesh(points, elements);

    % Assemble mass, stiffness, and T2-relaxation matrices in compartment
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
    K_cmpts{icmpt} = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));
    R_cmpts{icmpt} = 1 / relaxation(icmpt) * M_cmpts{icmpt};
    
    % Assemble moment matrices (coordinate weighted mass matrices)
    for idim = 1:3
        Jx_cmpts{idim}{icmpt} = mass_matrixP1_3D(elements', volumes', points(idim, :)');
    end
    
    % Create initial conditions (enforce complex values)
    rho_cmpts{icmpt} = complex(initial_density(icmpt)) * ones(npoint_cmpts(icmpt), 1);
end

% Create global mass, stiffness, relaxation, flux, and moment matrices (sparse)
disp("Coupling FEM matrices");
M = blkdiag(M_cmpts{:});
K = blkdiag(K_cmpts{:});
R = blkdiag(R_cmpts{:});
Jx = cellfun(@(J) blkdiag(J{:}), Jx_cmpts, "UniformOutput", false);
Q_blocks = assemble_flux_matrix(femesh.points, femesh.facets);
Q = couple_flux_matrix(femesh, setup.pde, Q_blocks, false);

% Global initial conditions
rho = vertcat(rho_cmpts{:});

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

    % Extract iteration inputs
    q = qvalues(iamp, iseq);
    b = bvalues(iamp, iseq);
    seq = sequences{iseq};
    g = directions(:, idir);

    % Get intervals based on the properties of the time profile
    [timelist, interval_str, timeprofile_str] = seq.intervals;

    % Number of intervals
    ninterval = length(timelist) - 1;

    % Assemble gradient direction dependent finite element matrix
    A = g(1) * Jx{1} + g(2) * Jx{2} + g(3) * Jx{3};

    % Base information about current iteration
    iteration_str = sprintf("Solving BTPDE of size %d using %s\n" ...
        + "  Direction %d of %d: g = [%.2f; %.2f; %.2f]\n" ...
        + "  Sequence  %d of %d: f = %s\n" ...
        + "  Amplitude %d of %d: q = %g, b = %g", ...
        sum(npoint_cmpts), solver_str, ...
        idir, ndirection, g, ...
        iseq, nsequence, seq, ...
        iamp, namplitude, q, b);

    % Initial conditions
    t = 0;
    y = rho;

    % Solve for each interval consecutively
    for iint = 1:ninterval

        % Display state of iterations
        fprintf("%s\n" ...
            + "  Interval %d of %d: I = %s, %s\n", ...
            iteration_str, ...
            iint, ninterval, interval_str(iint), timeprofile_str(iint));

        % Midpoint of interval
        interval_midpoint = (timelist(iint) + timelist(iint + 1)) / 2;

        % Right hand side Jacobian
        J = -(K + Q + R + 1i * seq.call(interval_midpoint) * q * A);

        % Factorize left hand side matrix
        [L, U, P, QQ, D] = lu(M - dt * theta * J);

        % Right hand side matrix
        E = M + dt * (1 - theta) * J;

        % Time step loop
        while t + dt < timelist(iint + 1)
            % Advance by dt
            % t
            y = QQ * (U \ (L \ (P * (D \ (E * y)))));
            t = t + dt;
        end

        % Adapt time step to advance to end of interval. This requires a new
        % LU decomposition
        dt_last = timelist(iint + 1) - t;
        E = M + dt_last * (1 - theta) * J;
        [L, U, P, QQ, D] = lu(M - dt_last * theta * J);
        y = QQ * (U \ (L \ (P * (D \ (E * y)))));
        t = t + dt_last;
    end

    % Split global solution into compartments
    mag = mat2cell(y, npoint_cmpts).';
    for icmpt = 1:ncompartment
        magnetization{icmpt, iall} = mag{icmpt};
    end
    signal(:, iall) = cellfun(@(M, y) sum(M * y, 1), M_cmpts, mag);

    % Store timing
    itertimes(iall) = toc(itertime);

end % gradient sequence iterations

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
