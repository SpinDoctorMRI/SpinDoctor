function results = solve_mf_cpu(femesh, setup, lap_eig, savepath, save_magnetization)
%SOLVE_MF_CPU Compute the solution to the BTPDE using Matrix Formalism on CPU.
%
%   SOLVE_MF_CPU(FEMESH, SETUP, LAP_EIG) solves the BTPDE using
%   Matrix Formalism and returns results.
%
%   SOLVE_MF_CPU(FEMESH, SETUP, LAP_EIG, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   SOLVE_MF_CPU(FEMESH, SETUP, LAP_EIG, SAVEPATH, SAVE_MAGNETIZATION) also omits saving
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
%   results: struct with fields. Split into the experiments for constant
%   direction vector sequences (const) and those with varying direction
%   (camino) from camino files. 
%   If only const or only camino sequences are present, then this
%   additional struct layer is removed.
%   
%       camino.magnetization: {ncompartment x nsequeunce}[npoint x 1]
%          Magnetization field at final timestep
%       camino.signal: [ncompartment x nsequence]
%           Compartmentwise total magnetization at final timestep
%       camino.signal_allcmpts: [nsequence x 1]
%           Total magnetization at final timestep
%       camino.itertimes: [nsequence x 1]
%           Computational time for each iteration
%       const.magnetization: {ncompartment x namplitude x nsequence x
%                       ndirection}[npoint x 1]
%           Magnetization field at final timestep
%       const.signal: [ncompartment x namplitude x nsequence x ndirection]
%           Compartmentwise total magnetization at final timestep
%       const.signal_allcmpts: [namplitude x nsequence x ndirection]
%           Total magnetization at final timestep
%       const.itertimes: [namplitude x nsequence x ndirection]
%           Computational time for each iteration
%       totaltime: [1 x 1]
%           Total computational time, including matrix assembly


% Measure time of function evaluation
starttime = tic;

% Check if a savepath has been provided (this triggers saving)
do_save = nargin >= nargin(@solve_mf_cpu) - 1;

% Define datatype, single precision helps reduce memory consumption
% and improve speed but degrade precision
if setup.mf.single
    dtype = @single;
else
    dtype = @double;
end

% Provide default value if not given
if nargin < nargin(@solve_mf_cpu)
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
const_sequences_ind = cellfun(@(x) ~isa(x,"SequenceCamino"),sequences,'UniformOutput',true);
nsequence_const = sum(const_sequences_ind);
sequences_const = sequences(const_sequences_ind);
const = struct;
const.magnetization = cell(ncompartment, namplitude, nsequence_const, ndirection);
const.signal = inf(ncompartment, namplitude, nsequence_const, ndirection);
const.signal_allcmpts = zeros(namplitude, nsequence_const, ndirection);
const.itertimes = zeros(namplitude, nsequence_const, ndirection);

nsequence_camino = sum(~const_sequences_ind);
camino = struct;
camino.magnetization = cell(ncompartment,nsequence_camino, 1);
camino.signal = inf(ncompartment, nsequence_camino);
camino.signal_allcmpts = zeros(nsequence_camino,1);
sequences_camino=sequences(~const_sequences_ind);
camino.itertimes = zeros(nsequence_camino, 1);
% Load if results are already available
if rerun && do_save
    fprintf("Load mf results from %s\n", savepath);
    totaltime_addition = 0;
    % Checking for const sequences
    for iseq = 1:nsequence_const
        % Extract iteration inputs
        seq = sequences_const{iseq};
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
                fprintf("Load mf for %s, %d/%d.\n", seq.string, iall, namplitude*ndirection);
                time_temp = totaltime_addition;
                try
                    savedata = mfile.(gradient_field);
                    const.signal(:, iamp, iseq, idir) = savedata.signal;
                    const.itertimes(iamp, iseq, idir) = savedata.itertimes;
                    totaltime_addition = totaltime_addition + savedata.itertimes;
                    if save_magnetization
                        const.magnetization(:, iamp, iseq, idir) = savedata.magnetization;
                    end
                catch
                    const.signal(:, iamp, iseq, idir) = inf;
                    const.itertimes(iamp, iseq, idir) = 0;
                    totaltime_addition = time_temp;
                    if save_magnetization
                        for icmpt = 1:ncompartment
                            const.magnetization{icmpt, iamp, iseq, idir} = [];
                        end
                    end
                    warning("mf: the saved data of experiment %s %s doesn't exist or is broken."...
                     + " Rerun simulation.", ...
                        seq.string, gradient_field);
                end
            end
        end
    end
    % Checking for camino sequences
    for iseq = 1:nsequence_camino
        seq = sequences_camino{iseq};
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        mfile = matfile(filename, "Writable", false);
        if hasfield(mfile,seq.string(true))
            % Load results
            fprintf("Load mf for %s \n", seq.string);
            time_temp = totaltime_addition;
            try
                savedata = mfile.(seq.string);
                camino.signal(:,iseq) = savedata.signal;
                camino.itertimes(iseq) = savedata.itertimes;
                totaltime_addition = totaltime_addition + savedata.itertimes;
                if save_magnetization
                    camino.magnetization(:, iseq) = savedata.magnetization;
                end
            catch e
                    disp(e)
                    camino.signal(:, iseq) = inf;
                    camino.itertimes( iseq) = 0;
                    totaltime_addition = time_temp;
                    if save_magnetization
                        for icmpt = 1:ncompartment
                            camnino.magnetization{icmpt, iseq} = [];
                        end
                    end
                    warning("mf: the saved data of experiment %s doesn't exist or is broken."...
                     + " Rerun simulation.", ...
                        seq.string);

            end
        end
    end
end

% Record unsaved experiments
no_result_flag_camino = permute(any(isinf(camino.signal), 1), [2 3 4 1]);
no_result_flag_const = permute(any(isinf(const.signal), 1), [2 3 4 1]);

if any(no_result_flag_camino, 'all') || any(no_result_flag_const, 'all')
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
            nu0 = dtype(funcs' * (M * rho));
            % Prepare LT2
            if no_relaxation
                LT2 = dtype(values');
            else
                LT2 = dtype(diag(values)+T2);
            end

            % save final laplace coefficient in nu_list
            nu_list_const = zeros(neig, namplitude, nsequence_const, ndirection, func2str(dtype));

            fprintf("Computing or loading MF magnetization for constant direction sequences for compartment %d " ...
                + "using %d eigenvalues.\n", ilapeig, neig);
            if any(no_result_flag_const, 'all')
                for idir = 1:ndirection
                    % Experiment parameters
                    ug = directions(:, idir);
                    
                    % Components of BT operator matrix
                    A = sum(moments .* shiftdim(ug, -2), 3);
                    A = dtype(A);
     
                    for iseq = 1:nsequence_const
                        % Experiment parameters
                        seq = sequences_const{iseq};
    
                        for iamp = 1:namplitude
                            % skip, if signal is already there
                            if all(~isinf(const.signal(:, iamp, iseq, idir)), 'all')
                                continue
                            else
                                no_result_flag_const(iamp, iseq, idir) = true;
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
                            iseq, nsequence_const, seq, ...
                            iamp, namplitude, g, q, b);
    
                            % Compute final laplace coefficient
                            nu = evolve_laplace_coef(nu0, seq, 1i*q*A, LT2, ninterval);
    
                            % Save final laplace coefficient in nu_list
                            nu_list_const(:, iamp, iseq, idir) = nu;
    
                            % Save computational time
                            const.itertimes(ilapeig, iamp, iseq, idir) = toc(itertime);
                        end
                    end
                end % iterations
    
                % Compute final magnetization
                nu_list_const = nu_list_const(:, no_result_flag_const);
                mag_const= funcs * double(nu_list_const);
                
                % Final magnetization coefficients in finite element nodal basis
                if save_magnetization
                    const.magnetization(ilapeig, no_result_flag_const) = ...
                        mat2cell(mag_const, size(mag_const, 1), ones(1, size(mag_const, 2)));
                end
                const.signal(ilapeig, no_result_flag_const) = sum(M * mag_const);
            end
            fprintf("Computing or loading MF magnetization for camino file sequences for compartment %d " ...
                + "using %d eigenvalues.\n", ilapeig, neig);
            if any(no_result_flag_camino, 'all')
                % save final laplace coefficient in nu_list
                nu_list = zeros(neig, nsequence_camino, func2str(dtype));
                for iseq = 1:nsequence_camino
                    % Experiment parameters
                    seq = sequences_camino{iseq};

                    % skip, if signal is already there
                    if all(~isinf(camino.signal(:, iseq, 1)), 'all')
                        continue
                    else
                        no_result_flag_camino( iseq, 1) = true;
                    end

                    % Measure iteration time
                    itertime = tic;
                    q = setup.gamma/1e6;
                    % Experiment parameters

                    % Run simulation if no result is saved or results are not available
                    % Display state of iterations
                    fprintf("Computing MF magnetization using %d eigenvalues\n" ...
                    + "  Sequence  %d of %d: f = %s\n", ...
                    neig,...
                    iseq, nsequence_camino, seq);

                    % Compute final laplace coefficient
                    nu = evolve_laplace_coef_direction_varying(nu0,seq, q,moments, LT2,dtype);
                    % Save final laplace coefficient in nu_list
                    nu_list(:, iseq) = nu;

                    % Save computational time
                    camino.itertimes(iseq) = toc(itertime);
                end

                % Compute final magnetization
                nu_list = nu_list(:, no_result_flag_camino);
                mag = funcs * double(nu_list);
                
                % Final magnetization coefficients in finite element nodal basis
                if save_magnetization
                    camino.magnetization(ilapeig, no_result_flag) = ...
                        mat2cell(mag, size(mag, 1), ones(1, size(mag, 2)));
                end
                camino.signal(ilapeig, no_result_flag) = sum(M * mag);
            end
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
        nu0 = dtype(funcs' * (M * rho));
        
        % Prepare LT2
        if no_relaxation
            LQT2 = dtype(LQ);
        else
            if size(LQ, 1) == 1
                LQT2 = dtype(diag(LQ)+T2);
            else
                LQT2 = dtype(LQ+T2);
            end
        end

        % save final laplace coefficient in nu_list
        nu_list_const = zeros(neig, namplitude, nsequence_const, ndirection, func2str(dtype));

        fprintf("Computing or loading MF magnetization using %d eigenvalues.\n", neig);
        if any(no_result_flag_const, 'all')
            for idir = 1:ndirection
                % Experiment parameters
                ug = directions(:, idir);
                
                % Components of BT operator matrix
                A = sum(moments .* shiftdim(ug, -2), 3);
                A = dtype(A);
    
                for iseq = 1:nsequence_const
                    % Experiment parameters
                    seq = sequences_const{iseq};
    
                   for iamp = 1:namplitude
                        % skip, if signal is already there
                        if all(~isinf(const.signal(:, iamp, iseq, idir)), 'all')
                            continue
                        else
                            no_result_flag_const(iamp, iseq, idir) = true;
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
                        iseq, nsequence_camino, seq, ...
                        iamp, namplitude, g, q, b);
    
                        % Compute final laplace coefficient
                        nu = evolve_laplace_coef(nu0, seq, 1i*q*A, LQT2, ninterval);
    
                        % Save final laplace coefficient in nu_list
                        nu_list_const(:, iamp, iseq, idir) = nu;
    
                        % Save computational time
                        const.itertimes(:, iamp, iseq, idir) = toc(itertime) * npoint_cmpts / sum(npoint_cmpts);
                    end
                end
            end % iterations
    
            % Compute final magnetization
            nu_list_const = nu_list_const(:, no_result_flag_const);
            mag_const= funcs * double(nu_list_const);
            
            % Final magnetization coefficients in finite element nodal basis
            idx = 1;
            allinds = [namplitude nsequence_const ndirection];
            for iall = 1:prod(allinds)
                % Extract indices
                [iamp, iseq, idir] = ind2sub(allinds, iall);
    
                if no_result_flag_const(iamp, iseq, idir)
                    mag_temp = mat2cell(mag_const(:, idx), npoint_cmpts);
                    if save_magnetization
                        const.magnetization(:, iamp, iseq, idir) = mag_temp;
                    end
                    const.signal(:, iamp, iseq, idir) = cellfun(@(M, m) sum(M * m), M_cmpts', mag_temp);
                    idx = idx + 1;
                end
            end
        end
    
        fprintf("Computing or loading MF magnetization for camino file sequences using %d eigenvalues.\n", neig);
        if any(no_result_flag_camino, 'all')
            % save final laplace coefficient in nu_list
            nu_list = zeros(neig, nsequence, func2str(dtype));
            for iseq = 1:nsequence_camino
                % Experiment parameters
                seq = sequences_camino{iseq};

                % skip, if signal is already there
                if all(~isinf(camino.signal(:, iseq, 1)), 'all')
                    continue
                else
                    no_result_flag_camino( iseq, 1) = true;
                end

                % Measure iteration time
                itertime = tic;
                q = setup.gamma/1e6;
                % Run simulation if no result is saved or results are not available
                % Display state of iterations
                fprintf("Computing MF magnetization using %d eigenvalues\n" ...
                + "  Sequence  %d of %d: f = %s\n",...
                neig, iseq, nsequence, seq);

                % Compute final laplace coefficient
                nu = evolve_laplace_coef_direction_varying(nu0,seq, q,moments, LQT2,dtype);
                % Save final laplace coefficient in nu_list
                nu_list(:, iseq) = nu;

                % Save computational time
                camino.itertimes(iseq) = toc(itertime) * npoint_cmpts / sum(npoint_cmpts);
            end
                % end
            % Compute final magnetization
            nu_list = nu_list(:, no_result_flag_camino);
            mag = funcs * double(nu_list);
            
            % Final magnetization coefficients in finite element nodal basis
            idx = 1;
            
            for iseq = 1:nsequence_camino
                % Extract indices
                if no_result_flag_camino(iseq, 1)
                    mag_temp = mat2cell(mag(:, idx), npoint_cmpts);
                    if save_magnetization
                        camino.magnetization(:,  iseq) = mag_temp;
                    end
                    camino.signal(:, iseq, 1) = cellfun(@(M, m) sum(M * m), M_cmpts', mag_temp);
                    idx = idx + 1;
                end
            end
        end
    end
end

if do_save && any(no_result_flag_const, 'all')
    for iseq = 1:nsequence_const
        seq = sequences_const{iseq};
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        fprintf("Save %s\n", filename);
        mfile = matfile(filename, "Writable", true);
        for iamp = 1:namplitude
            for idir = 1:ndirection
                if no_result_flag_const(iamp, iseq, idir)
                    % Extract iteration inputs
                    data = struct;
                    data.q = qvalues(iamp,iseq);
                    data.b = bvalues(iamp, iseq);
                    data.ug = directions(:, idir);
                    data.g = gvalues(iamp,iseq);
                    data.signal = const.signal(:,iamp, iseq, idir);
                    data.itertimes = const.itertimes(iamp,iseq,idir);
                    if save_magnetization
                        data.magnetization = const.magnetization(:,iamp,iseq,idir);
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

if do_save && any(no_result_flag_camino, 'all')
    for iseq = 1:nsequence_camino
        seq = sequences_camino{iseq};
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        fprintf("Save %s\n", filename);
        mfile = matfile(filename, "Writable", true);
       
        if no_result_flag_camino(iseq, 1)
            % Extract iteration inputs
            data = struct;
            data.seq= seq;
            data.signal = camino.signal(:,  iseq);
            data.sequence = seq;
            data.itertimes = camino.itertimes(iseq);
            if save_magnetization
                data.magnetization = camino.magnetization(:, iseq);
            end

            % Save results to MAT-file
            gradient_field = sprintf("%s",seq);
            mfile.(gradient_field) = data;
        end
    end
end
% Total magnetization (sum over compartments)
camino.signal_allcmpts(:) = sum(camino.signal, 1);
const.signal_allcmpts(:) = sum(const.signal, 1);
% Create output structure
results.totaltime = toc(starttime);
if nsequence_camino == 0
    results.signal_weighted = const.signal./femesh.volumes;
    results.signal_allcmpts_weighted = const.signal_allcmpts;
    results.signal = const.signal;
    results.signal_allcmpts = const.signal_allcmpts;

    results.itertimes = const.itertimes;
    if save_magnetization
        results.magnetization = const.magnetization;
        results.magnetization_avg = average_magnetization(const.magnetization);
    end
elseif nsequence_const == 0
    results.signal = camino.signal;
    results.signal_allcmpts = camino.signal_allcmpts;
    results.signal_weighted = camino.signal./femesh.volumes;
    results.signal_allcmpts_weighted = camino.signal_allcmpts/femesh.total_volume;
    results.itertimes = camino.itertimes;
    if save_magnetization
        results.magnetization = camino.magnetization;
        results.magnetization_avg = average_magnetization(camino.magnetization);
    end

else
    results.camino.signal = camino.signal;
    results.camino.signal_allcmpts = camino.signal_allcmpts;
    
    results.camino.signal_weighted = camino.signal./femesh.volumes;
    results.camino.signal_allcmpts_weighted = camino.signal_allcmpts/femesh.total_volume;
    results.camino.itertimes = camino.itertimes;

    results.const.signal = const.signal;
    results.const.signal_allcmpts = const.signal_allcmpts;
    results.const.signal = const.signal;
    results.const.signal_allcmpts = const.signal_allcmpts;
    results.const.itertimes = const.itertimes;

    if save_magnetization
        results.camino.magnetization = camino.magnetization;
        results.camino.magnetization_avg = average_magnetization(camino.magnetization);
        results.const.magnetization = const.magnetization;
        results.const.magnetization_avg = average_magnetization(const.magnetization);
    end

end 
% Display function evaluation time
toc(starttime);
