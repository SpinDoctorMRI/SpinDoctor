function results = solve_mf_cpu_no_lap(setup, savepath)
%SOLVE_MF_CPU_NO_LAP Compute the solution to the BTPDE using Matrix Formalism on CPU.
%
%
%   SOLVE_MF_CPU_NO_LAP(SETUP, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/<DMRIINFO>/<MF_INFO>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   setup: struct
%   savepath: path string
%
%   results: struct with fields. Split into the experiments for constant
%   direction vector sequences (const) and those with varying direction
%   (camino) from camino files. 
%   If only const or only camino sequences are present, then this
%   additional struct layer is removed.
%   
%       camino.signal: [ncompartment x nsequence]
%           Compartmentwise total magnetization at final timestep
%       camino.signal_allcmpts: [nsequence x 1]
%           Total magnetization at final timestep
%       camino.itertimes: [nsequence x 1]
%           Computational time for each iteration
%       const.signal: [ncompartment x namplitude x nsequence x ndirection]
%           Compartmentwise total magnetization at final timestep
%       const.signal_allcmpts: [namplitude x nsequence x ndirection]
%           Total magnetization at final timestep
%       const.itertimes: [namplitude x nsequence x ndirection]
%           Computational time for each iteration
%       totaltime: [1 x 1]
%           Total computational time, including matrix assembly

disp("Running solve_mf_cpu_no_lap, work in progress.");

% Measure time of function evaluation
starttime = tic;

% Check if a savepath has been provided (this triggers saving)
do_save = nargin >= nargin(@solve_mf_cpu_no_lap) - 1;

% Define datatype, single precision helps reduce memory consumption
% and improve speed but degrade precision
if setup.mf.single
    dtype = @single;
else
    dtype = @double;
end

if isfield(setup.mf,'rerun')
    rerun = setup.mf.rerun;
else
    rerun = false;
end
% Extract domain parameters
relaxation = setup.pde.relaxation;
no_relaxation = all(isinf(relaxation));
multi_lap_eig = setup.ncompartment > 1; %length(lap_eig) > 1;
if multi_lap_eig
    error("Not currently supported for multi_lap_eig");
end


save_magnetization=false;

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
% FEM_save_path = make_FEM_save_path(setup.pde,setup.mf,multi_lap_eig,savepath);

if do_save
    % Folder for saving and loading is the same
    savepath = add_mf_str_savepath(savepath,setup.mf );
    loadpath = savepath;
else
    % No folder for saving, create folder for loading.
    savepath = ""; loadpath = create_savepath(setup,"mf_no_lap_eig");
    loadpath = add_mf_str_savepath(loadpath,setup.mf );
end

if ~isfolder(loadpath)
    error("Finite element matrices must have already been saved to %s\nRun solve_mf with laplacian eigenfunctions.",savepath);
end
FEM_save_path = sprintf("%s/FEM_lap_eig.mat",loadpath);

% Initialize output arguments
const_sequences_ind = cellfun(@(x) ~isa(x,"SequenceCamino"),sequences,'UniformOutput',true);
nsequence_const = sum(const_sequences_ind);
sequences_const = sequences(const_sequences_ind);
const = struct;
const.signal = inf(ncompartment, namplitude, nsequence_const, ndirection);
const.signal_allcmpts = zeros(namplitude, nsequence_const, ndirection);
const.itertimes = zeros(namplitude, nsequence_const, ndirection);

nsequence_camino = sum(~const_sequences_ind);
camino = struct;
camino.signal = inf(ncompartment, nsequence_camino);
camino.signal_allcmpts = zeros(nsequence_camino,1);
sequences_camino=sequences(~const_sequences_ind);
camino.itertimes = zeros(nsequence_camino, 1);

% Load if results are already available
totaltime_addition = 0;
if ~rerun && do_save
    fprintf("Load mf results from %s\n", savepath);
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
                    
                catch
                    const.signal(:, iamp, iseq, idir) = inf;
                    const.itertimes(iamp, iseq, idir) = 0;
                    totaltime_addition = time_temp;
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
            catch e
                disp(e)
                camino.signal(:, iseq) = inf;
                camino.itertimes( iseq) = 0;
                totaltime_addition = time_temp;
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
                
                const.signal(ilapeig, no_result_flag_const) = sum(M * mag_const);
            end

            if any(no_result_flag_camino, 'all')
                fprintf("Computing MF magnetization for camino file sequences for compartment %d " ...
                + "using %d eigenvalues.\n", ilapeig, neig);
                % save final laplace coefficient in nu_list
                nu_list_camino = zeros(neig, nsequence_camino, func2str(dtype));
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
                    nu_list_camino(:, iseq) = nu;

                    % Save computational time
                    camino.itertimes(iseq) = toc(itertime);
                end

                % Compute final magnetization
                nu_list_camino = nu_list_camino(:, no_result_flag_camino);
                mag = funcs * double(nu_list_camino);
                
                camino.signal(ilapeig, no_result_flag) = sum(M * mag);
            end
        end % lapeig iterations
    else
        % Prepare mass, density, moments and relaxation matrices
        FEM = load(FEM_save_path);
        moments =FEM.moments;
        % npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

        % Coefficients of initial spin density in Laplace eigenfunction basis
        nu0 = FEM.nu0;
        % Prepare LT2
        LQT2 = FEM.LQT2;
        neig = size(LQT2,2);
        nu2signal = FEM.nu2signal;
        clear FEM;
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
                        iseq, nsequence_const, seq, ...
                        iamp, namplitude, g, q, b);
    
                        % Compute final laplace coefficient
                        nu = evolve_laplace_coef(nu0, seq, 1i*q*A, LQT2, ninterval);
    
                        % Save final laplace coefficient in nu_list
                        nu_list_const(:, iamp, iseq, idir) = nu;
    
                        % Save computational time
                        const.itertimes(:, iamp, iseq, idir) = toc(itertime);% * npoint_cmpts / sum(npoint_cmpts);
                    end
                end
            end % iterations
    
            % Compute final magnetization
            nu_list_const = nu_list_const(:, no_result_flag_const);            
            % Final magnetization coefficients in finite element nodal basis
            allinds = [namplitude nsequence_const ndirection];
            idx= 1;
            for iall = 1:prod(allinds)
                % Extract indices
                [iamp, iseq, idir] = ind2sub(allinds, iall);
    
                if no_result_flag_const(iamp, iseq, idir)
                    const.signal(:, iamp, iseq, idir) = nu2signal*nu_list_const(:, idx) ;
                    idx = idx +1;
                end
            end
        end
    
        fprintf("Computing or loading MF magnetization for camino file sequences using %d eigenvalues.\n", neig);
        if any(no_result_flag_camino, 'all')
            % save final laplace coefficient in nu_list
            nu_list_camino = zeros(neig, nsequence, func2str(dtype));
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
                nu_list_camino(:, iseq) = nu;

                % Save computational time
                camino.itertimes(iseq) = toc(itertime);% * npoint_cmpts / sum(npoint_cmpts);
            end
                % end
            % Compute final magnetization
            nu_list_camino = nu_list_camino(:, no_result_flag_camino);
            % mag = funcs * double(nu_list_camino);
            
            % Final magnetization coefficients in finite element nodal basis
            idx = 1;
            for iseq = 1:nsequence_camino
                % Extract indices
                if no_result_flag_camino(iseq, 1)

                    camino.signal(:, iseq, 1) = nu2signal*nu_list_camino(:,idx);
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
totaltime = totaltime_addition + toc(starttime);
results = merge_results(camino,const,nsequence_camino,nsequence_const,totaltime,save_magnetization);

% Display function evaluation time
toc(starttime);
