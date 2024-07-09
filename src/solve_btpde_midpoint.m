function results = solve_btpde_midpoint(femesh, setup, savepath, save_magnetization)
%SOLVE_BTPDE_MIDPOINT Solve the Bloch-Torrey partial differential equation.
%   This solver uses a theta timestepping method (generalized midpoint), and
%   only works for interval-wise constant time profiles (PGSE, DoublePGSE).
%
%   SOLVE_BTPDE_MIDPOINT(FEMESH, SETUP) solves the BTPDE and returns results.
%
%   SOLVE_BTPDE_MIDPOINT(FEMESH, SETUP, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/btpde_midpoint/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads the results
%   instead of solving for that iteration.
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
%       camino.magnetization: {ncompartment x nsequence}[npoint x 1]
%          Magnetization field at final timestep
%       const.signal: [ncompartment x namplitude x nsequence x ndirection]
%           Compartmentwise total magnetization at final timestep
%       const.signal_allcmpts: [namplitude x nsequence x ndirection]
%           Total magnetization at final timestep
%       const.itertimes: [namplitude x nsequence x ndirection]
%           Computational time for each iteration
%       const.magnetization: {ncompartment x namplitude x nsequence x
%                       ndirection}[npoint x 1]
%           Magnetization field at final timestep
%       totaltime: [1 x 1]
%           Total computational time, including matrix assembly

% TEMPORARY. Camino file sequences not yet implemented for this solver.
% const_ind = cellfun(@(x) ~isa(x,"SequenceCamino"),setup.gradient.sequences,'UniformOutput',true);
% if ~all(const_ind,'all')
%     warning("Currently %s does not support camino file sequences. \n Solving only for non-camino sequences",mfilename);
%     setup.gradient.sequences = setup.gradient.sequences(const_ind);
%     setup.nsequence = sum(const_ind);
% end

% Measure function evaluation time
starttime = tic;

% Check if a save path has been provided (this toggers saving)
do_save = nargin >= nargin(@solve_btpde_midpoint) - 1;

% Provide default value
if nargin < nargin(@solve_btpde_midpoint)
    save_magnetization = true;
end
if isfield(setup.btpde_midpoint, 'rerun')
    rerun = setup.btpde_midpoint.rerun;
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
theta = setup.btpde_midpoint.implicitness;
dt = setup.btpde_midpoint.timestep;
solver_str = sprintf("theta rule (theta = %g, dt = %g)", theta, dt);

% Check that sequences are compatible
% assert(all(cellfun(@(f) isa(f, "PGSE") || isa(f, "DoublePGSE"), sequences)), ...
%     "Only constant sequences are currently supported.");

% Sizes
ncompartment = setup.ncompartment;
namplitude = setup.namplitude;
% nsequence = setup.nsequence;
ndirection = setup.ndirection;

% Number of points in each compartment
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

if do_save
    % Folder for saving
    savepath = sprintf( ...
        "%s/theta%g_dt%g", ...
        savepath, theta, dt ...
    );
    if ~isfolder(savepath)
        mkdir(savepath)
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
totaltime_addition = 0;

% Check if results are already available
if ~rerun && do_save
    fprintf("Load btpde results from %s\n", savepath);
    
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
                fprintf("Load btpde for %s, %d/%d.\n", seq.string, iall, namplitude*ndirection);
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
                    warning("btpde: the saved data of experiment %s %s doesn't exist or is broken."...
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
            fprintf("Load btpde for %s \n", seq.string);
            time_temp = totaltime_addition;
            try
                savedata = mfile.(seq.string);
                camino.signal(:,iseq) = savedata.signal;
                camino.itertimes(iseq) = savedata.itertimes;
                totaltime_addition = totaltime_addition + savedata.itertimes;
                    if save_magnetization
                        camino.magnetization(:, iseq) = savedata.magnetization;
                    end
            catch
                    camino.signal(:, iseq) = inf;
                    const.itertimes( iseq) = 0;
                    totaltime_addition = time_temp;
                    if save_magnetization
                        for icmpt = 1:ncompartment
                            camino.magnetization{icmpt, iseq} = [];
                        end
                    end
                    warning("btpde: the saved data of experiment %s doesn't exist or is broken."...
                     + " Rerun simulation.", ...
                        seq.string);

            end
        end
    end
end




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


% Iterate over gradient amplitudes, sequences and directions. If the Matlab
% PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
% parallel, otherwise it should work like a normal loop. If that is not the
% case, replace the `parfor` keyword by the normal `for` keyword.

if any(isinf(const.signal), 'all') 
    disp("Simulating for constant direction vector sequences")
    % Cartesian indices (for parallel looping with linear indices)
    allinds_const = [namplitude nsequence_const ndirection];

    % Iterate over gradient amplitudes, sequences and directions. If the Matlab
    % PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
    % parallel, otherwise it should work like a normal loop. If that is not the
    % case, replace the `parfor` keyword by the normal `for` keyword.

    % Temporarily save results in temp_store to avoid I/O error
    temp_store = cell(allinds_const);

    % Check if Parallel Computing Toolbox is licensed
    if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
        parpool('local', [1, 2048]);
    end
    signal = const.signal;
    magnetization = const.magnetization;
    itertimes = const.itertimes;
    parfor iall = 1:prod(allinds_const)
        % skip, if signal is already there
        if all(~isinf(signal(:, iall)), 'all')
            continue
        end
        
        % Measure iteration time
        itertime = tic;
        % Extract Cartesian indices
        [iamp, iseq, idir] = ind2sub(allinds_const, iall);
        % Extract iteration inputs
        seq = sequences_const{iseq};
        q = qvalues(iamp, iseq);
        b = bvalues(iamp, iseq);
        ug = directions(:, idir);
        g = gvalues(iamp, iseq);

        % Solve the pde for a pulse sequence.
        mag = btpde_midpoint_seq_const(seq,ug,q,b,g,rho,dt, ...
            theta,M,Jx,K,Q,R,idir,iseq,iamp,ndirection, ...
            nsequence_const,namplitude,npoint_cmpts,solver_str);
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
    const.signal= signal;
    const.sequences = sequences;
    const.magnetization = magnetization;
    const.itertimes = itertimes;
    % Save data into
    if do_save
        for iseq = 1:nsequence_const
            seq = sequences_const{iseq};
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
                        % const.itertimes(iamp,iseq,idir) = temp_store{iamp, iseq, idir}.itertimes;
                    end
                end
            end
        end

    end
end

if any(isinf(camino.signal),'all')    
    disp("Simulating for camino file sequences")
    q = setup.gamma/1e6;
    signal = camino.signal;
    sequences = sequences_camino;
    magnetization = camino.magnetization;
    itertimes = camino.itertimes;
    parfor iseq = 1:nsequence_camino
        % skip, if signal is already there
        if all(~isinf(signal(:, iseq)), 'all')
            continue
        end
        
        % Measure iteration time
        itertime = tic;
    
        % Extract Cartesian indices
        % Extract iteration inputs
        seq = sequences{iseq};
        mag = btpde_midpoint_seq_camino(seq,q,rho,dt,theta,M,Jx,K,Q,R,iseq,nsequence_camino,npoint_cmpts,solver_str);
        % Split global solution into compartments
        mag = mat2cell(mag, npoint_cmpts).';
        signal(:, iseq) = cellfun(@(M, y) sum(M * y, 1), M_cmpts, mag);
    
        % Store timing
        itertimes(iseq) = toc(itertime);
        
        if do_save
            data = struct;
            data.signal = signal(:, iseq);
            data.itertimes = itertimes(iseq);
            if save_magnetization
                data.magnetization = mag;
            end
    
            % Save iteration results
            temp_store{iseq} = data;
        end
        
        % Store magnetization
        if save_magnetization
            magnetization(:, iseq) = mag;
        end
    end % iterations
    camino.signal= signal;
    camino.sequences = sequences;
    camino.magnetization = magnetization;
    camino.itertimes = itertimes;

    % Saving data
    if do_save
        for iseq = 1:nsequence_camino
            seq = sequences_camino{iseq};
            filename = sprintf("%s/%s.mat", savepath, seq.string(true));
            fprintf("Save %s\n", filename);
            mfile = matfile(filename, "Writable", true);
            if ~isempty(temp_store{iseq})
                mfile.(seq.string()) = temp_store{iseq};
            end
        end

    end

end
% Total magnetization (sum over compartments)
camino.signal_allcmpts(:) = sum(camino.signal, 1);
const.signal_allcmpts(:) = sum(const.signal, 1);

totaltime = totaltime_addition + toc(starttime);
results = merge_results(camino,const,nsequence_camino,nsequence_const,totaltime,save_magnetization);

% Display function evaluation time
toc(starttime);
