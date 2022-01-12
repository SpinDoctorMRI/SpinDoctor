function results = solve_mf(femesh, setup, lap_eig, savepath, save_magnetization)
%SOLVE_MF Compute the solution to the BTPDE using Matrix Formalism.
%
%   SOLVE_MF(FEMESH, SETUP, LAP_EIG) solves the BTPDE using
%   Matrix Formalism and returns results.
%
%   SOLVE_MF(FEMESH, SETUP, LAP_EIG, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   SOLVE_MF(FEMESH, SETUP, LAP_EIG, SAVEPATH, SAVE_MAGNETIZATION) also omits saving
%   or loading the magnetization field if SAVE_MAGNETIZATION is set to FALSE.
%
%   femesh: struct
%   setup: struct
%   lap_eig: struct with fields
%       values: double(neig, 1)
%       funcs: double(npoint, neig)
%   savepath (optional): path string
%   save_magnetization (optinal): logical. Defaults to true.
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


% Measure time of function evaluation
starttime = tic;

% Check if a savepath has been provided (this toggers saving)
do_save = nargin >= nargin(@solve_mf) - 1;

% Provide default value if not given
if nargin < nargin(@solve_btpde)
    save_magnetization = true;
end
rerun = setup.mf.rerun;

% Extract domain parameters
initial_density = setup.pde.initial_density;
relaxation = setup.pde.relaxation;
no_relaxation = all(isinf(relaxation));
multi_lap_eig = length(lap_eig) > 1;

% Extract experiment parameters
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
gvalues = setup.gradient.gvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;
ninterval = setup.mf.ninterval;

% Sizes
ncompartment = setup.ncompartment;
namplitude = setup.namplitude;
nsequence = setup.nsequence;
ndirection = setup.ndirection;

if do_save
    % Folder for saving
    mf_str = sprintf("neig_max%g_lengthscale_min%.4f_ninterval%d", ...
        setup.mf.neig_max, setup.mf.length_scale, setup.mf.ninterval);
    if ~isinf(setup.mf.neig_max)
        % if neig_max is inf, mf.eigs doesn't exist or is removed.
        mf_str = mf_str + sprintf("_md5_%s", DataHash(setup.mf.eigs, 10));
    end
    savepath = fullfile(savepath, mf_str);
    if ~isfolder(savepath)
        mkdir(savepath);
    end
else
    savepath = "";
end

% Initialize output arguments
magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(ncompartment, namplitude, nsequence, ndirection);

% Assemble mass matrix in each compartment (for spatial integration)
disp("Setting up FEM matrices");
M_cmpts = cell(1, ncompartment);
R_cmpts = cell(1, ncompartment);
Jx_cmpts = repmat({cell(1, ncompartment)}, 1, 3);
rho_cmpts = cell(1, ncompartment);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [~, volumes] = get_volume_mesh(points, elements);

    % Assemble mass and T2-relaxation matrices in compartment
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
    if no_relaxation
        R_cmpts{icmpt} = 0;
    else
        R_cmpts{icmpt} = 1 / relaxation(icmpt) * M_cmpts{icmpt};
    end
    % Assemble moment matrices (coordinate weighted mass matrices)
    for idim = 1:3
        Jx_cmpts{idim}{icmpt} = mass_matrixP1_3D(elements', volumes', points(idim, :)');
    end
    % Create initial conditions (enforce complex values)
    rho_cmpts{icmpt} = complex(initial_density(icmpt)) * ones(size(femesh.points{icmpt}, 2), 1);
end

% Cartesian indices (for looping with linear indices)
ainds = [namplitude ndirection];
% Record unsaved experiments
no_result_flag = false(namplitude, nsequence, ndirection);

if multi_lap_eig
    for ilapeig = 1:length(lap_eig)
        % Extract eigenvalues, eigenfunctions
        L = diag(lap_eig(ilapeig).values);
        funcs = lap_eig(ilapeig).funcs;
        neig = length(lap_eig(ilapeig).values);

        % Prepare mass, density, moments and relaxation matrices
        M = M_cmpts{ilapeig};
        rho = rho_cmpts{ilapeig};
        % Compute first order moments of eigenfunction products
        Jx = cellfun(@(J) J{ilapeig}, Jx_cmpts, "UniformOutput", false);
        moments = zeros(neig, neig, 3);
        for idim = 1:3
            moments(:, :, idim) = funcs' * Jx{idim} * funcs;
        end
        % Compute T2-weighted Laplace mass matrix
        if no_relaxation
            T2 = 0;
        else
            R = R_cmpts{ilapeig};
            T2 = funcs' * R * funcs;
        end
        % Coefficients of initial spin density in Laplace eigenfunction basis
        nu0 = funcs' * (M * rho);

        fprintf("Computing or loading MF magnetization for compartment %d " ...
            + "using %d eigenvalues.\n", ilapeig, neig);

        for iseq = 1:nsequence
            % Sequence parameters
            seq = sequences{iseq};
            % File for saving or loading iteration results
            filename = sprintf("%s/%s.mat", savepath, seq.string(true));
            mfile = matfile(filename, "Writable", false);

            for ia = 1:prod(ainds)
                % Measure iteration time
                itertime = tic;
    
                % Extract indices
                [iamp, idir] = ind2sub(ainds, ia);
    
                % Experiment parameters
                q = qvalues(iamp, iseq);
                b = bvalues(iamp, iseq);
                ug = directions(:, idir);
                g = gvalues(iamp, iseq);
    
                gradient_field = gradient_fieldstring(ug, b);
                no_result = true;
    
                % Check if results are already available
                if ~rerun && do_save && hasfield(mfile, gradient_field)
                    % Load results
                    try
                        data = mfile.(gradient_field);
                        signal(ilapeig, iamp, iseq, idir) = data.signal(ilapeig);
                        if save_magnetization
                            magnetization{ilapeig, iamp, iseq, idir} = data.magnetization{ilapeig};
                        end
                        itertimes(ilapeig, iamp, iseq, idir) = data.itertimes(ilapeig);
                        no_result = false;
                    catch
                        no_result = true;
                        warning("mf: the saved data of experiment %s %s is broken. Rerun simulation.", ...
                            seq.string, gradient_field);
                    end
                end
                no_result_flag(ia) = no_result_flag(ia) | no_result;
    
                % Run simulation if no result is saved or results are not available
                if no_result
                    % Display state of iterations
                    fprintf("Computing MF magnetization using %d eigenvalues\n" ...
                    + "  Direction %d of %d: ug = [%.2f; %.2f; %.2f]\n" ...
                    + "  Sequence  %d of %d: f = %s\n" ...
                    + "  Amplitude %d of %d: g = %g, q = %g, b = %g\n", ...
                    neig, ...
                    idir, ndirection, ug, ...
                    iseq, nsequence, seq, ...
                    iamp, namplitude, g, q, b);
    
                    % Components of BT operator matrix
                    A = sum(moments .* shiftdim(ug, -2), 3);
    
                    % Compute final magnetization (in Laplace basis)
                    nu = evolve_laplace_coef(nu0, seq, 1i * q * A, L + T2, ninterval);
    
                    % Final magnetization coefficients in finite element nodal basis
                    mag = funcs * nu;
                    if save_magnetization
                        magnetization{ilapeig, iamp, iseq, idir} = mag;
                    end
                    signal(ilapeig, iamp, iseq, idir) = sum(M * mag);
                    itertimes(ilapeig, iamp, iseq, idir) = toc(itertime);
                end
            end % iterations
        end % sequences
    end % lapeig iterations
else
    % Extract eigenvalues, eigenfunctions
    L = diag(lap_eig.values);
    funcs = lap_eig.funcs;
    neig = length(lap_eig.values);
    % Number of points in each compartment
    npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

    % Prepare mass, density, moments and relaxation matrices
    M = blkdiag(M_cmpts{:});
    rho = vertcat(rho_cmpts{:});
    % Compute first order moments of eigenfunction products
    Jx = cellfun(@(J) blkdiag(J{:}), Jx_cmpts, "UniformOutput", false);
    moments = zeros(neig, neig, 3);
    for idim = 1:3
        moments(:, :, idim) = funcs' * Jx{idim} * funcs;
    end
    % Compute T2-weighted Laplace mass matrix
    if no_relaxation
        T2 = 0;
    else
        R = blkdiag(R_cmpts{:});
        T2 = funcs' * R * funcs;
    end
    % Coefficients of initial spin density in Laplace eigenfunction basis
    nu0 = funcs' * (M * rho);

    fprintf("Computing or loading MF magnetization using %d eigenvalues.\n", neig);
    for iseq = 1:nsequence
        % Sequence parameters
        seq = sequences{iseq};
        % File name for saving or loading iteration results
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        mfile = matfile(filename, "Writable", false);

        for ia = 1:prod(ainds)
            % Measure iteration time
            itertime = tic;
    
            % Extract indices
            [iamp, idir] = ind2sub(ainds, ia);
    
            % Experiment parameters
            q = qvalues(iamp, iseq);
            b = bvalues(iamp, iseq);
            ug = directions(:, idir);
            g = gvalues(iamp, iseq);
    
            gradient_field = gradient_fieldstring(ug, b);
            no_result = true;
    
            % Check if results are already available
            if ~rerun && do_save && hasfield(mfile, gradient_field)
                % Load results
                fprintf("Load mf for %s, %d/%d.\n", seq.string, ia, namplitude*nsequence*ndirection);
                try
                    data = mfile.(gradient_field);
                    signal(:, iamp, iseq, idir) = data.signal;
                    if save_magnetization
                        magnetization(:, iamp, iseq, idir) = data.magnetization;
                    end
                    itertimes(:, iamp, iseq, idir) = data.itertimes;
                    no_result = false;
                catch
                    no_result = true;
                    warning("mf: the saved data of experiment %s %s is broken. Rerun simulation.", ...
                        seq.string, gradient_field);
                end
            end
            no_result_flag(ia) = no_result_flag(ia) | no_result;
    
            % Run simulation if no result is saved or results are not available
            if no_result
                % Display state of iterations
                fprintf("Computing MF magnetization using %d eigenvalues\n" ...
                + "  Direction %d of %d: ug = [%.2f; %.2f; %.2f]\n" ...
                + "  Sequence  %d of %d: f = %s\n" ...
                + "  Amplitude %d of %d: g = %g, q = %g, b = %g\n", ...
                neig, ...
                idir, ndirection, ug, ...
                iseq, nsequence, seq, ...
                iamp, namplitude, g, q, b);
    
                % Components of BT operator matrix
                A = sum(moments .* shiftdim(ug, -2), 3);
    
                % Compute final magnetization (in Laplace basis)
                nu = evolve_laplace_coef(nu0, seq, 1i * q * A, L + T2, ninterval);
    
                % Final magnetization coefficients in finite element nodal basis
                mag = funcs * nu;
                mag = mat2cell(mag, npoint_cmpts);
                if save_magnetization
                    magnetization(:, iamp, iseq, idir) = mag;
                end
                signal(:, iamp, iseq, idir) = cellfun(@(M, m) sum(M * m), M_cmpts', mag);
                itertimes(:, iamp, iseq, idir) = toc(itertime) * npoint_cmpts / sum(npoint_cmpts);
            end
        end % iterations
    end % sequences
end

if do_save && any(no_result_flag, 'all')
    for iseq = 1:nsequence
        seq = sequences{iseq};
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        fprintf("Save %s\n", filename);
        mfile = matfile(filename, "Writable", true);
        for iamp = 1:namplitude
            for idir = 1:ndirection
                if no_result_flag(iamp, iseq, idir)
                    % Extract iteration inputs
                    data = struct;
                    data.q = qvalues(iamp, iseq);
                    data.b = bvalues(iamp, iseq);
                    data.ug = directions(:, idir);
                    data.g = gvalues(iamp, iseq);
                    data.signal = signal(:, iamp, iseq, idir);
                    data.itertimes = itertimes(:, iamp, iseq, idir);
                    if save_magnetization
                        data.magnetization = magnetization(:, iamp, iseq, idir);
                    end

                    % Save results to MAT-file
                    gradient_field = gradient_fieldstring(data.ug, data.b);
                    mfile.(gradient_field) = data;

                    % dMRI signal is centrosymmetric
                    data.ug = -data.ug;
                    % convert negative zeros to positive zeros
                    data.ug(data.ug == 0) = +0;
                    gradient_field = gradient_fieldstring(data.ug, data.b);
                    if ~hasfield(mfile, gradient_field)
                        mfile.(gradient_field) = data;
                    end
                end
            end
        end
    end
end

% Compute total signal (sum over compartments)
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime);
if save_magnetization
    results.magnetization = magnetization;
    results.magnetization_avg = average_magnetization(magnetization);
end

% Display function evaluation time
toc(starttime);
