function results = solve_mf_gpu(femesh, setup, lap_eig, savepath, save_magnetization)
%SOLVE_MF_GPU Compute the solution to the BTPDE using Matrix Formalism on GPU.
%
%   SOLVE_MF_GPU(FEMESH, SETUP, LAP_EIG) solves the BTPDE using
%   Matrix Formalism and returns results.
%
%   SOLVE_MF_GPU(FEMESH, SETUP, LAP_EIG, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   SOLVE_MF_GPU(FEMESH, SETUP, LAP_EIG, SAVEPATH, SAVE_MAGNETIZATION) also omits saving
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

% TEMPORARY. Camino file sequences not yet implemented for this solver.
const_ind = cellfun(@(x) ~isa(x,"SequenceCamino"),setup.gradient.sequences,'UniformOutput',true);
if ~all(const_ind,'all')
    warning("Currently %s does not support camino file sequences. \n Solving only for non-camino sequences",mfilename);
    setup.gradient.sequences = setup.gradient.sequences(const_ind);
    setup.nsequence = sum(const_ind);
end

% Measure time of function evaluation
starttime = tic;

% Check if a savepath has been provided (this triggers saving)
do_save = nargin >= nargin(@solve_mf_gpu) - 1;

% Define datatype, single precision helps reduce memory consumption
% and improve speed but degrade precision
if setup.mf.single
    dtype = @single;
else
    dtype = @double;
end

% Provide default value if not given
if nargin < nargin(@solve_mf_gpu)
    save_magnetization = true;
end
rerun = setup.mf.rerun;

% Extract domain parameters
initial_density = setup.pde.initial_density;
relaxation = setup.pde.relaxation;
no_relaxation = all(isinf(relaxation));
zero_permeability = all(setup.pde.permeability==0);
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
    mf_str = sprintf("neig%g_ls%.4f", ...
        setup.mf.neig_max, setup.mf.length_scale);
    if setup.mf.surf_relaxation
        mf_str = "surf_relaxation_" + mf_str;
    end
    if setup.mf.single
        mf_str = mf_str + "_single";
    end
    if ~isinf(setup.mf.neig_max)
        % if neig_max is inf, mf.eigs doesn't exist or is removed.
        mf_str = mf_str + sprintf("_%s", DataHash(setup.mf.eigs, 6));
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
signal = inf(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(ncompartment, namplitude, nsequence, ndirection);

% Load if results are already available
if ~rerun && do_save
    fprintf("Load mf results from %s\n", savepath);

    % Iterate over gradient amplitudes, sequences and directions
    for iseq = 1:nsequence
        seq = sequences{iseq};
        % Load results
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        mfile = matfile(filename, "Writable", false);

        for ia = 1:prod([namplitude, ndirection])
            % Extract Cartesian indices
            [iamp, idir] = ind2sub([namplitude, ndirection], ia);

            % Extract iteration inputs
            b = bvalues(iamp, iseq);
            ug = directions(:, idir);

            % File name for saving or loading iteration results
            gradient_field = gradient_fieldstring(ug, b);

            % Check if results are already available
            if hasfield(mfile, gradient_field)
                % Load results
                fprintf("Load mf for %s, %d/%d.\n", ...
                    seq.string, ia, namplitude*ndirection);

                try
                    data = mfile.(gradient_fieldstring(ug, b));
                    signal(:, iamp, iseq, idir) = data.signal;
                    itertimes(:, iamp, iseq, idir) = data.itertimes;
                    if save_magnetization
                        magnetization(:, iamp, iseq, idir) = data.magnetization;
                    end
                catch
                    signal(:, iamp, iseq, idir) = inf;
                    itertimes(:, iamp, iseq, idir) = 0;
                    if save_magnetization
                        for icmpt = 1:ncompartment
                            magnetization{icmpt, iamp, iseq, idir} = [];
                        end
                    end
                    warning("mf: the saved data of experiment %s %s doesn't exist or is broken."...
                     + " Rerun simulation.", ...
                        seq.string, gradient_field);
                end
            end
        end
    end
end

% Record unsaved experiments
no_result_flag = permute(any(isinf(signal), 1), [2 3 4 1]);

if any(no_result_flag, 'all')
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

    if setup.mf.surf_relaxation && ~zero_permeability
        % construct and assemble flux matrix
        Q_blocks = assemble_flux_matrix(femesh.points, femesh.facets);
        Q = couple_flux_matrix(femesh, setup.pde, Q_blocks, false);
    end

    % Select GPU device
    if isnumeric(setup.mf.gpu)
        gpuDevice(setup.mf.gpu);
    else
        gpuDevice(); 
    end

    if multi_lap_eig && ~setup.mf.surf_relaxation
        for ilapeig = 1:length(lap_eig)
            % Extract eigenvalues, eigenfunctions
            values = lap_eig(ilapeig).values;
            funcs = lap_eig(ilapeig).funcs;
            neig = length(values);

            % Prepare mass, density, moments and T2-relaxation matrices
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
            nu0 = gpuArray(dtype(funcs' * (M * rho)));
            % Prepare LT2
            if no_relaxation
                LT2 = gpuArray(dtype(values'));
            else
                LT2 = gpuArray(dtype(diag(values)+T2));
            end

            % save final laplace coefficient in nu_list
            nu_list = zeros(neig, namplitude, nsequence, ndirection, func2str(dtype));

            fprintf("Computing or loading MF magnetization for compartment %d " ...
                + "using %d eigenvalues.\n", ilapeig, neig);
            for idir = 1:ndirection
                % Experiment parameters
                ug = directions(:, idir);
                
                % Components of BT operator matrix
                A = sum(moments .* shiftdim(ug, -2), 3);
                A = gpuArray(dtype(A));
 
                for iseq = 1:nsequence
                    % Experiment parameters
                    seq = sequences{iseq};

                    for iamp = 1:namplitude
                        % skip, if signal is already there
                        if all(~isinf(signal(:, iamp, iseq, idir)), 'all')
                            continue
                        else
                            no_result_flag(iamp, iseq, idir) = true;
                        end

                        % Measure iteration time
                        itertime = tic;

                        % Experiment parameters
                        q = qvalues(iamp, iseq);
                        b = bvalues(iamp, iseq);
                        g = gvalues(iamp, iseq);

                        % Run simulation if no result is saved or results are not available
                        % Display state of iterations
                        fprintf("Computing MF magnetization using %d eigenvalues\n" ...
                        + "  Direction %d of %d: ug = [%.2f; %.2f; %.2f]\n" ...
                        + "  Sequence  %d of %d: f = %s\n" ...
                        + "  Amplitude %d of %d: g = %g, q = %g, b = %g\n", ...
                        neig, ...
                        idir, ndirection, ug, ...
                        iseq, nsequence, seq, ...
                        iamp, namplitude, g, q, b);

                        % Compute final laplace coefficient
                        nu = evolve_laplace_coef(nu0, seq, 1i*q*A, LT2, ninterval);

                        % Save final laplace coefficient in nu_list
                        nu_list(:, iamp, iseq, idir) = gather(nu);

                        % Save computational time
                        itertimes(ilapeig, iamp, iseq, idir) = toc(itertime);
                    end
                end
            end % iterations

            % release GPU memory
            gpuDevice([]);

            % Compute final magnetization
            nu_list = nu_list(:, no_result_flag);
            mag = compute_mag_gpu(dtype(funcs), dtype(nu_list), setup.mf.gpu);
            
            % Final magnetization coefficients in finite element nodal basis
            if save_magnetization
                magnetization(ilapeig, no_result_flag) = ...
                    mat2cell(mag, size(mag, 1), ones(1, size(mag, 2)));
            end
            signal(ilapeig, no_result_flag) = sum(M * mag);
        end % lapeig iterations
    else
        % Prepare mass, density, moments and relaxation matrices
        M = blkdiag(M_cmpts{:});
        rho = vertcat(rho_cmpts{:});

        % Number of points in each compartment
        npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

        % Extract eigenvalues, eigenfunctions
        if setup.mf.surf_relaxation
            disp("Matrix formalism: surface relaxation is on.");
            if multi_lap_eig
                % Assemble eigenvalues and eigenfunctions of all decoupled compartments
                [values, funcs] = assemble_lapeig(lap_eig, M, npoint_cmpts);
            else
                values = lap_eig.values;
                funcs = lap_eig.funcs;
            end
            % Surface relaxation matrix
            if zero_permeability
                LQ = values';
            else
                surface_relaxation = funcs' * Q * funcs;
                LQ = diag(values) + surface_relaxation;
            end
        else
            LQ = lap_eig.values';
            funcs = lap_eig.funcs;
        end
        neig = size(LQ, 2);

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
        nu0 = gpuArray(dtype(funcs' * (M * rho)));
        % Prepare LT2
        if no_relaxation
            LQT2 = gpuArray(dtype(LQ));
        else
            if size(LQ, 1) == 1
                LQT2 = gpuArray(dtype(diag(LQ)+T2));
            else
                LQT2 = gpuArray(dtype(LQ+T2));
            end
        end

        % save final laplace coefficient in nu_list
        nu_list = zeros(neig, namplitude, nsequence, ndirection, func2str(dtype));

        fprintf("Computing or loading MF magnetization using %d eigenvalues.\n", neig);
        for idir = 1:ndirection
            % Experiment parameters
            ug = directions(:, idir);
            
            % Components of BT operator matrix
            A = sum(moments .* shiftdim(ug, -2), 3);
            A = gpuArray(dtype(A));

            for iseq = 1:nsequence
                % Experiment parameters
                seq = sequences{iseq};

                for iamp = 1:namplitude   
                    % skip, if signal is already there
                    if all(~isinf(signal(:, iamp, iseq, idir)), 'all')
                        continue
                    else
                        no_result_flag(iamp, iseq, idir) = true;
                    end

                    % Measure iteration time
                    itertime = tic;

                    % Experiment parameters
                    q = qvalues(iamp, iseq);
                    b = bvalues(iamp, iseq);
                    g = gvalues(iamp, iseq);

                    % Run simulation if no result is saved or results are not available
                    % Display state of iterations
                    fprintf("Computing MF magnetization using %d eigenvalues\n" ...
                    + "  Direction %d of %d: ug = [%.2f; %.2f; %.2f]\n" ...
                    + "  Sequence  %d of %d: f = %s\n" ...
                    + "  Amplitude %d of %d: g = %g, q = %g, b = %g\n", ...
                    neig, ...
                    idir, ndirection, ug, ...
                    iseq, nsequence, seq, ...
                    iamp, namplitude, g, q, b);

                    % Compute final laplace coefficient
                    nu = evolve_laplace_coef(nu0, seq, 1i*q*A, LQT2, ninterval);

                    % Save final laplace coefficient in nu_list
                    nu_list(:, iamp, iseq, idir) = gather(nu);

                    % Save computational time
                    itertimes(:, iamp, iseq, idir) = toc(itertime) * npoint_cmpts / sum(npoint_cmpts);
                end
            end
        end % iterations

        % release GPU memory
        gpuDevice([]);

        % Compute final magnetization
        nu_list = nu_list(:, no_result_flag);
        mag = compute_mag_gpu(dtype(funcs), dtype(nu_list), setup.mf.gpu);

        % Final magnetization coefficients in finite element nodal basis
        idx = 1;
        allinds = [namplitude nsequence ndirection];
        for iall = 1:prod(allinds)
            % Extract indices
            [iamp, iseq, idir] = ind2sub(allinds, iall);

            if no_result_flag(iamp, iseq, idir)
                mag_temp = mat2cell(mag(:, idx), npoint_cmpts);
                if save_magnetization
                    magnetization(:, iamp, iseq, idir) = mag_temp;
                end
                signal(:, iamp, iseq, idir) = cellfun(@(M, m) sum(M * m), M_cmpts', mag_temp);
                idx = idx + 1;
            end
        end
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
end
