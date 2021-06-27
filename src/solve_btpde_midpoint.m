function results = solve_btpde_midpoint(femesh, setup, savepath, save_magnetization)
%SOLVE_BTPDE_MIDPOINT Solve the Bloch-Torrey partial differential equation.
%   This solver uses a theta timestepping method (generalized midpoint), and
%   only works for interval-wise constant time profiles (PGSE, DoublePGSE).
%
%   SOLVE_BTPDE_MIDPOINT(FEMESH, SETUP) solves the BTPDE and returns results.
%
%   SOLVE_BTPDE_MIDPOINT(FEMESH, SETUP, SAVEPATH) saves the results of each
%   iteration at "<SAVEPATH>/BTPDE_<SOLVEROPTIONS>/<ITERATIONINFO>.MAT". If an
%   iteration file is already present, the solver loads the results instead of
%   solving for that iteration.
%
%   SOLVE_BTPDE_MIDPOINT(FEMESH, SETUP, SAVEPATH, SAVE_MAGNETIZATION) also omits
%   saving or loading the magnetization field if SAVE_MAGNETIZATION is set to
%   FALSE.
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
%   femesh: struct
%   setup: struct
%   savepath (optional): string
%   save_magnetization (optional): logical. Defaults to true.
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

% Check if a save path has been provided (this toggers saving)
do_save = nargin >= nargin(@solve_btpde_midpoint) - 1;

% Provide default value
if nargin < nargin(@solve_btpde_midpoint)
    save_magnetization = true;
end

% Extract domain parameters
diffusivity = setup.pde.diffusivity;
relaxation = setup.pde.relaxation;
initial_density = setup.pde.initial_density;

% Extract experiment parameters
values = setup.gradient.values;
amptype = setup.gradient.values_type;
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

if do_save
    % Folder for saving
    savepath = sprintf( ...
        "%s/btpde_midpoint_theta%g_dt%g_magnetization%d", ...
        savepath, theta, dt, save_magnetization ...
    );
    if ~isfolder(savepath)
        mkdir(savepath)
    end
else
    savepath = "";
end

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);
totaltime_addition = 0;

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
    amp = values(iamp);
    q = qvalues(iamp, iseq);
    b = bvalues(iamp, iseq);
    seq = sequences{iseq};
    g = directions(:, idir);
    
    % File name for saving or loading iteration results
    filename = sprintf("%s/%s.mat", savepath, gradient_string(amp, amptype, seq, g));
    
    % Check if results are already available
    if do_save && isfile(filename)
        % Load results
        fprintf("Load %s\n", filename);
        mfile = matfile(filename, "Writable", false);
        signal(:, iall) = mfile.signal;
        itertimes(iall) = mfile.itertime;
        totaltime_addition = totaltime_addition + mfile.itertime;
        if save_magnetization
            mag = mfile.magnetization;
        end
    else
        % Get intervals based on the properties of the time profile
        [timelist, interval_str, timeprofile_str] = seq.intervals;

        % Number of intervals
        ninterval = length(timelist) - 1;

        % Assemble gradient direction dependent finite element matrix
        A = g(1) * Jx{1} + g(2) * Jx{2} + g(3) * Jx{3};

        % Initial conditions
        t = 0;
        y = rho;

        % Solve for each interval consecutively
        for iint = 1:ninterval
        
            % Display state of iterations
            fprintf( ...
                join([
                    "Solving BTPDE of size %d using %s"
                    "  Direction %d of %d: g = [%.2f; %.2f; %.2f]"
                    "  Sequence  %d of %d: f = %s"
                    "  Amplitude %d of %d: q = %g, b = %g"
                    "  Interval  %d of %d: I = %s, %s\n"
                ], newline), ...
                sum(npoint_cmpts), solver_str, ...
                idir, ndirection, g, ...
                iseq, nsequence, seq, ...
                iamp, namplitude, q, b, ...
                iint, ninterval, interval_str(iint), timeprofile_str(iint) ...
            );
            
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
        signal(:, iall) = cellfun(@(M, y) sum(M * y, 1), M_cmpts, mag);

        % Store timing
        itertimes(iall) = toc(itertime);
        
        if do_save
            % Save iteration results
            fprintf("Save %s\n", filename);
            mfile = matfile(filename, "Writable", true);
            mfile.signal = signal(:, iall);
            mfile.itertime = itertimes(iall);
            if save_magnetization
                mfile.magnetization = mag;
            end
        end
        
    end % load or save variables
    
    % Store magnetization
    if save_magnetization
        for icmpt = 1:ncompartment
            magnetization{icmpt, iall} = mag{icmpt};
        end
    end
end % gradient sequence iterations

% Total magnetization (sum over compartments)
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.magnetization = magnetization;
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime) + totaltime_addition;

% Display function evaluation time
toc(starttime);
