function results = solve_hadc(femesh, setup, savepath)
%SOLVE_HADC Compute the ADC from HADC model.
%
%   SOLVE_HADC(FEMESH, SETUP) solves the HADC and returns results.
%
%   SOLVE_HADC(FEMESH, SETUP, SAVEPATH) saves the results of each iteration at
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/hadc/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%   If a result is already present in the iteration file, the solver loads
%   the results instead of solving for that iteration.
%
%   femesh: struct
%   setup: struct
%   savepath (optional): string
%
%   results: struct with fields
%       adc: double(ncompartment, nsequence, ndirection)
%       adc_allcmpts: double(nsequence, ndirection)
%       itertimes: double(ncompartment, nsequence, ndirection)
%       totaltime: double(ncompartment, nsequence, ndirection)


% Check if a save path has been provided (this toggers saving)
do_save = nargin == nargin(@solve_hadc);
if isfield(setup.hadc, 'rerun')
    rerun = setup.hadc.rerun;
else
    rerun = false;
end

% Measure function evaluation time
starttime = tic;

% Extract experiment parameters
directions = setup.gradient.directions;
sequences = setup.gradient.sequences;
reltol = setup.hadc.reltol;
abstol = setup.hadc.abstol;
solve_ode = setup.hadc.ode_solver;
solver_str = func2str(solve_ode);

% Extract domain parameters
diffusivity = setup.pde.diffusivity;
initial_density = setup.pde.initial_density;
compartments = setup.pde.compartments;

% Sizes
ncompartment = length(compartments);
nsequence = length(sequences);
ndirection = size(directions, 2);

% Number of points in each compartment
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

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
adc = zeros(ncompartment, nsequence, ndirection);
adc_allcmpts = zeros(nsequence, ndirection);
itertimes = zeros(ncompartment, nsequence, ndirection);
totaltime_addition = 0;

% Assemble finite element matrices
M = cell(1, ncompartment);
K = cell(1, ncompartment);
G = cell(1, ncompartment);
rho = cell(1, ncompartment);
volumes = zeros(1, ncompartment);
for icmpt = 1:ncompartment
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [volumes(icmpt), fevolumes] = get_volume_mesh(points, elements);

    % Assemble finite element matrices
    M{icmpt} = mass_matrixP1_3D(elements', fevolumes');
    K{icmpt} = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));

    % Compute surface integrals
    G{icmpt} = zeros(npoint_cmpts(icmpt), 3);
    for iboundary = 1:femesh.nboundary
        facets = femesh.facets{icmpt, iboundary};
        if ~isempty(facets)
            % Get facet normals
            [~, ~, normals] = get_surfacenormal_mesh(points, elements, facets);

            % Surface normal weighted flux matrix (in each canonical direction)
            for dim = 1:3
                Q = flux_matrixP1_3D(facets', points', normals(dim, :)');
                G{icmpt}(:, dim) = G{icmpt}(:, dim) + sum(Q, 2);
            end
        end
    end

    % Initial conditions
    rho{icmpt} = zeros(npoint_cmpts(icmpt), 1);
end

% Cartesian indices (for parallel looping with linear indices)
allinds = [ncompartment nsequence ndirection];

% Iterate over compartments and gradient sequences and directions. If the Matlab
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
    % Measure iteration time
    itertime = tic;

    % Extract indices
    [icmpt, iseq, idir] = ind2sub(allinds, iall);

    % Extract parameters for iteration
    seq = sequences{iseq};
    ug = directions(:, idir);

    % File name for saving or loading iteration results
    filename = sprintf("%s/%s.mat", savepath, seq.string(true));
    mfile = matfile(filename, "Writable", false);
    gradient_field = sprintf("cmpt%d_", icmpt) + gradient_fieldstring(ug);
    no_result = true;
    
    % Check if results are already available
    if ~rerun && do_save && hasfield(mfile, gradient_field)
        % Load results
        fprintf("Load hadc %d/%d.\n", iall, prod(allinds));
        try
            data = mfile.(gradient_field);
            adc(iall) = data.adc;
            itertimes(iall) = data.itertimes;
            totaltime_addition = totaltime_addition + data.itertimes;
            no_result = false;
        catch
            no_result = true;
            warning("hadc: the saved data of experiment %s %s is broken. Rerun simulation.", ...
                seq.string, gradient_field);
        end
    end

    % Run simulation if no result is saved or results are not available
    if no_result
        % Free diffusivity in direction ug
        D0 = ug' * diffusivity(:, :, icmpt) * ug;

        % Compute surface integrals in gradient direction
        surfint = G{icmpt} * (diffusivity(:, :, icmpt) * ug);

        % Display state of iterations
        fprintf( ...
            join([
                "Solving HADC model of size %d using %s:"
                "  Direction %d of %d: ug = [%.2f; %.2f; %.2f]"
                "  Sequence  %d of %d: f = %s"
                "  Compartment %d of %d: %s\n"
            ], newline), ...
            sum(npoint_cmpts), solver_str, ...
            idir, ndirection, ug, ...
            iseq, nsequence, seq, ...
            icmpt, ncompartment, compartments(icmpt) ...
        );
        
        % We only keep two points in tlist to output solution at all time
        % steps
        tlist = [0 seq.echotime];

        % Create ODE function and Jacobian from matrices
        [ode_function, Jacobian] = hadc_functions(K{icmpt}, surfint, seq);

        % Set parameters for ODE solver
        options = odeset( ...
            "Mass", M{icmpt}, ...
            "AbsTol", abstol, ...
            "RelTol", reltol, ...
            "Vectorized", "on", ...
            "Stats", "off", ...
            "Jacobian", Jacobian ...
        );

        % Solve ODE, keep all time steps (for integral)
        [tvec, y] = solve_ode(ode_function, tlist, rho{icmpt}, options);
        tvec = tvec';
        y = y';

        % Integral over compartment boundary
        hvec = surfint' * y / volumes(icmpt);

        % HADC (free diffusivity minus correction)
        a = trapz(tvec, seq.integral(tvec) .* hvec) / seq.bvalue_no_q;
        adc(iall) = D0 - a;

        % Computational time
        itertimes(iall) = toc(itertime);
        
        if do_save
            data.ug = ug;
            data.adc = adc(iall);
            data.itertimes = itertimes(iall);
            
            % Save iteration results
            temp_store{iall} = data;
        end
    end
end

if do_save
    for iseq = 1:nsequence
        seq = sequences{iseq};
        filename = sprintf("%s/%s.mat", savepath, seq.string(true));
        fprintf("Save %s\n", filename);
        mfile = matfile(filename, "Writable", true);
        for icmpt = 1:ncompartment
            for idir = 1:ndirection
                if ~isempty(temp_store{icmpt, iseq, idir})
                    % Extract iteration inputs
                    ug = directions(:, idir);

                    % Save results to MAT-file
                    gradient_field = sprintf("cmpt%d_", icmpt) + gradient_fieldstring(ug);
                    mfile.(gradient_field) = temp_store{icmpt, iseq, idir};

                    % adc is centrosymmetric
                    ug = -ug;
                    % convert negative zeros to positive zeros
                    ug(ug == 0) = +0;
                    gradient_field = sprintf("cmpt%d_", icmpt) + gradient_fieldstring(ug);
                    if ~hasfield(mfile, gradient_field)
                        temp_store{icmpt, iseq, idir}.ug = ug;
                        mfile.(gradient_field) = temp_store{icmpt, iseq, idir};
                    end
                end
            end
        end
    end
end

% Compute total HADC (weighted sum over compartments)
weights = initial_density .* volumes;
weights = weights / sum(weights);
adc_allcmpts(:) = sum(weights' .* adc, 1);

% Create output structure
results.adc = adc;
results.adc_allcmpts = adc_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime) + totaltime_addition;

% Display function evaluation time
toc(starttime);
