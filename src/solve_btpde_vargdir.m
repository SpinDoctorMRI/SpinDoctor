function results = solve_btpde_vargdir(femesh,setup,savepath,save_magnetization)
%SOLVE_BTPDE Solve the Bloch-Torrey partial differential equation.
%
%   SOLVE_BTPDE(FEMESH, SETUP) solves the BTPDE and returns results.
%
%   SOLVE_BTPDE(FEMESH, SETUP, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/BTPDE_<SOLVEROPTIONS>/<ITERATIONINFO>.MAT". If an
%   iteration file is already present, the solver loads the results instead of
%   solving for that iteration.
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
reltol = setup.btpde.reltol;
abstol = setup.btpde.abstol;
solve_ode = setup.btpde.ode_solver;
solver_str = func2str(solve_ode);

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
        "%s/btpde_abstol%g_reltol%g_magnetization%d", ...
        savepath, abstol, reltol, save_magnetization ...
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

for iall = 1:prod(allinds)
    
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
        
        %%%%%%%%%% For CoreyBaron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Overriding the timelist
        timelist = setup.timelist_vargdir;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        % Number of intervals
        ninterval = length(timelist) - 1;
        
        
        % Initial magnetization
        mag = rho;
        
        % Solve for each interval consecutively
        for iint = 1:ninterval
          
            %%%%%%%%%% For CoreyBaron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%% Commented this line out that uses a fixed g
            % Assemble gradient direction dependent finite element matrix          
            %J = g(1) * Jx{1} + g(2) * Jx{2} + g(3) * Jx{3};
            
            %%% Define vargdir
            g_vargdir = setup.vargdir(iint,:);
            %%% Use vargdir here           
            J = g_vargdir(1) * Jx{1} + g_vargdir(2) * Jx{2} + g_vargdir(3) * Jx{3};
            
            %%% Define these quantities below for the new vargdir sequence, need them for
            %%% printing later in the code
            interval_str(iint) = num2str(iint);           
            timeprofile_str(iint) = timeprofile_str(1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Add a third point to the interval, so that the ODE solver does not
            % store the magnetization for all time steps during the solve. If
            % there were only two points in the interval, the ODE solver would
            % store all time steps. This would require a lot of memory,
            % especially during parfor iterations
            interval_midpoint = (timelist(iint) + timelist(iint + 1)) / 2;
            time_list_interval = [timelist(iint), interval_midpoint, timelist(iint + 1)];
            
            %%%%%%%%%% For CoreyBaron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Modified display to show the correct g at each time step
            %%% q is between 0 and 1 (0 is for simulating b=0, q=1 is for simulating 
            %%% the uploaded sequence)
            % Display state of iterations
            fprintf( ...
                join([
                "Solving BTPDE of size %d using %s:"
                "  Direction %d of %d: g = [%.2g; %.2g; %.2g]"
                "  Sequence  %d of %d: f = %s"
                "  Amplitude %d of %d: q = %g, b = %g"
                "  Interval  %d of %d: I = %s, %s\n"
                ], newline), ...
                sum(npoint_cmpts), solver_str, ...
                idir, ndirection, g_vargdir, ...
                iseq, nsequence, seq, ...
                iamp, namplitude, q, b, ...
                iint, ninterval, interval_str(iint), timeprofile_str(iint) ...
                );
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%% For CoreyBaron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Modified btpde_functions_interval function so that the
            %%% Jacobian is constant for this choice, the new function is
            %%% called btpde_functions_interval_vargdir
            % Create new ODE functions on given interval
            [ode_function, Jacobian] = btpde_functions_interval_vargdir( ...
                K, Q, R, J, q, seq, interval_midpoint);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
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
end % iterations

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
