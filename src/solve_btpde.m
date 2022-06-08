function results = solve_btpde(femesh, setup, savepath, save_magnetization)
%SOLVE_BTPDE Solve the Bloch-Torrey partial differential equation.
%
%   SOLVE_BTPDE(FEMESH, SETUP) solves the BTPDE and returns results.
%
%   SOLVE_BTPDE(FEMESH, SETUP, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   SOLVE_BTPDE(FEMESH, SETUP, SAVEPATH, SAVE_MAGNETIZATION) also omits saving
%   or loading the magnetization field if SAVE_MAGNETIZATION is set to FALSE.
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
do_save = nargin >= nargin(@solve_btpde) - 1;

% Provide default value
if nargin < nargin(@solve_btpde)
    save_magnetization = true;
end
if isfield(setup.btpde, 'rerun')
    rerun = setup.btpde.rerun;
else
    rerun = false;
end

% Extract domain parameters
diffusivity = setup.pde.diffusivity;
relaxation = setup.pde.relaxation;
initial_density = setup.pde.initial_density;

% Extract experiment parameters
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
gvalues = setup.gradient.gvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;
reltol = setup.btpde.reltol;
abstol = setup.btpde.abstol;
solve_ode = setup.btpde.ode_solver;
solver_str = func2str(solve_ode);

% Sizes
ncompartment = setup.ncompartment;
namplitude = setup.namplitude;
nsequence = setup.nsequence;
ndirection = setup.ndirection;

if do_save
    % Folder for saving
    savepath = sprintf( ...
        "%s/%s_abstol%g_reltol%g", ...
        savepath, solver_str, abstol, reltol ...
    );
    if ~isfolder(savepath)
        mkdir(savepath)
    end
else
    savepath = "";
end

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = inf(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);
totaltime_addition = 0;

% Check if results are already available
if ~rerun && do_save
    fprintf("Load btpde results from %s\n", savepath);

    for iseq = 1:nsequence
        % Extract iteration inputs
        seq = sequences{iseq};
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        mfile = matfile(filename, "Writable", false);

        for iall = 1:prod([namplitude, ndirection])   
            % Extract Cartesian indices
            [iamp, idir] = ind2sub([namplitude, ndirection], iall);
        
            % Extract iteration inputs
            b = bvalues(iamp, iseq);
            ug = directions(:, idir);
        
            % File name for saving or loading iteration results
            gradient_field = gradient_fieldstring(ug, b);
        
            % Check if results are already available
            if hasfield(mfile, gradient_field)
                % Load results
                fprintf("Load btpde for %s, %d/%d.\n", seq.string, iall, namplitude*ndirection);
                time_temp = totaltime_addition;
                try
                    savedata = mfile.(gradient_field);
                    signal(:, iamp, iseq, idir) = savedata.signal;
                    itertimes(iamp, iseq, idir) = savedata.itertimes;
                    totaltime_addition = totaltime_addition + savedata.itertimes;
                    if save_magnetization
                        magnetization(:, iamp, iseq, idir) = savedata.magnetization;
                    end
                catch
                    signal(:, iamp, iseq, idir) = inf;
                    itertimes(iamp, iseq, idir) = 0;
                    totaltime_addition = time_temp;
                    if save_magnetization
                        for icmpt = 1:ncompartment
                            magnetization{icmpt, iamp, iseq, idir} = [];
                        end
                    end
                    warning("btpde: the saved data of experiment %s %s doesn't exist or is broken."...
                     + " Rerun simulation.", ...
                        seq.string, gradient_field);
                end
            end
        end
    end
end

if any(isinf(signal), 'all')
    % Number of points in each compartment
    npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

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

    % Set parameters for ODE solver
    options_template = odeset( ...
        "Mass", M, ...
        "AbsTol", abstol, ...
        "RelTol", reltol, ...
        "Vectorized", "on", ...
        "Stats", "off", ...
        "MassSingular", "no" ...
    );

    % Cartesian indices (for parallel looping with linear indices)
    allinds = [namplitude nsequence ndirection];

    % Iterate over gradient amplitudes, sequences and directions. If the Matlab
    % PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
    % parallel, otherwise it should work like a normal loop. If that is not the
    % case, replace the `parfor` keyword by the normal `for` keyword.

    % Temporarily save results in temp_store to avoid I/O error
    temp_store = cell(allinds);

    % Check if Parallel Computing Toolbox is licensed
    if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
        parpool('local', [1, 2048]);
    end

    parfor iall = 1:prod(allinds)
        % skip, if signal is already there
        if all(~isinf(signal(:, iall)), 'all')
            continue
        end
        
        % Measure iteration time
        itertime = tic;

        % Extract Cartesian indices
        [iamp, iseq, idir] = ind2sub(allinds, iall);

        % Extract iteration inputs
        seq = sequences{iseq};
        q = qvalues(iamp, iseq);
        b = bvalues(iamp, iseq);
        ug = directions(:, idir);
        g = gvalues(iamp, iseq);

        % Get intervals based on the properties of the time profile
        [timelist, interval_str, timeprofile_str] = seq.intervals;

        % Number of intervals
        ninterval = length(timelist) - 1;

        % Assemble gradient direction dependent finite element matrix
        J = ug(1) * Jx{1} + ug(2) * Jx{2} + ug(3) * Jx{3};

        % Initial magnetization
        mag = rho;

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
            fprintf( ...
                join([
                    "Solving BTPDE of size %d using %s:"
                    "  Direction %d of %d: ug = [%.2f; %.2f; %.2f]"
                    "  Sequence  %d of %d: f = %s"
                    "  Amplitude %d of %d: g = %g, q = %g, b = %g"
                    "  Interval  %d of %d: I = %s, %s\n"
                ], newline), ...
                sum(npoint_cmpts), solver_str, ...
                idir, ndirection, ug, ...
                iseq, nsequence, seq, ...
                iamp, namplitude, g, q, b, ...
                iint, ninterval, interval_str(iint), timeprofile_str(iint) ...
            );

            % Create new ODE functions on given interval
            [ode_function, Jacobian] = btpde_functions_interval( ...
                K, Q, R, J, q, seq, interval_midpoint);

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
        signal(:, iall) = cellfun(@(M, y) sum(M * y, 1), M_cmpts, mag);

        % Store timing
        itertimes(iall) = toc(itertime);
        
        if do_save
            data = struct;
            data.b = b;
            data.q = q;
            data.g = g;
            data.ug = ug;
            data.signal = signal(:, iall);
            data.itertimes = itertimes(iall);
            if save_magnetization
                data.magnetization = mag;
            end

            % Save iteration results
            temp_store{iall} = data;
        end
        
        % Store magnetization
        if save_magnetization
            magnetization(:, iall) = mag;
        end
    end % iterations

    if do_save
        for iseq = 1:nsequence
            seq = sequences{iseq};
            filename = sprintf("%s/%s.mat", savepath, seq.string(true));
            fprintf("Save %s\n", filename);
            mfile = matfile(filename, "Writable", true);
            for iamp = 1:namplitude
                for idir = 1:ndirection
                    if ~isempty(temp_store{iamp, iseq, idir})
                        % Extract iteration inputs
                        b = bvalues(iamp, iseq);
                        ug = directions(:, idir);

                        % Save results to MAT-file
                        gradient_field = gradient_fieldstring(ug, b);
                        mfile.(gradient_field) = temp_store{iamp, iseq, idir};

                        % dMRI signal is centrosymmetric
                        ug = -ug;
                        % convert negative zeros to positive zeros
                        ug(ug == 0) = +0;
                        gradient_field = gradient_fieldstring(ug, b);
                        if ~hasfield(mfile, gradient_field)
                            temp_store{iamp, iseq, idir}.ug = ug;
                            mfile.(gradient_field) = temp_store{iamp, iseq, idir};
                        end
                    end
                end
            end
        end
    end
end

% Total magnetization (sum over compartments)
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime) + totaltime_addition;
if save_magnetization
    results.magnetization = magnetization;
    results.magnetization_avg = average_magnetization(magnetization);
end

% Display function evaluation time
toc(starttime);
