function results = solve_hadc(femesh, params_domain, experiment)
%SOLVE_HADC Compute the ADC from HADC model.
%
% 	femesh: struct
% 	params_domain: struct
% 	experiment: struct
%
%   results: struct with fields
%       adc: double(ncompartment, nsequence, ndirection)
%      	adc_allcmpts: double(nsequence, ndirection)
%      	ctime: double(ncompartment, nsequence, ndirection)


% Measure function evaluation time
starttime = tic;

% Extract experiment parameters
sequences = experiment.sequences;
reltol = experiment.hadc.reltol;
abstol = experiment.hadc.abstol;
solve_ode = experiment.hadc.ode_solver;
solver_str = func2str(solve_ode);

% Extract domain parameters
diffusivity = params_domain.diffusivity;
initial_density = params_domain.initial_density;
compartments = params_domain.compartments;

% Take direction averaged diffusivity from tensor (trace)
diffusivity_scalar = (diffusivity(1, 1, :) + diffusivity(2, 2, :) + diffusivity(3, 3, :)) / 3;
diffusivity_scalar = shiftdim(diffusivity_scalar, 1);

% Extract HARDI points
dir_points = experiment.directions.points;
dir_inds = experiment.directions.indices;
opposite = experiment.directions.opposite;

% Sizes
ncompartment = length(compartments);
nsequence = length(sequences);
ndirection = experiment.ndirection;
ndirection_unique = length(dir_inds);

% Number of points in each compartment
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);

% Initialize output arguments
adc = zeros(ncompartment, nsequence, ndirection);
adc_allcmpts = zeros(nsequence, ndirection);
itertimes = zeros(ncompartment, nsequence, ndirection);

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
    for dim = 1:3
        G{icmpt} = zeros(npoint_cmpts(icmpt), 3);
    end
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
parfor iall = 1:prod(allinds)

    % Measure iteration time
    itertime = tic;

    % Extract indices
    [icmpt, iseq, idir] = ind2sub(allinds, iall);

    % Do not bother solving if opposite direction has already been computed
    if idir <= ndirection_unique

        % Extract parameters for iteration
        seq = sequences{iseq};
        g = dir_points(:, idir);
        
        % Free diffusivity in direction g
        D0 = g' * diffusivity(:, :, icmpt) * g;

        % Compute surface integrals in gradient direction
        surfint = G{icmpt} * (diffusivity(:, :, icmpt) * g);

        % Display state of iterations
        fprintf("Solving HADC equation using %s\n" ...
            + "  Direction   %d of %d: g = [%.2f; %.2f; %.2f]\n" ...
            + "  Sequence    %d of %d: f = %s\n" ...
            + "  Compartment %d of %d: %s\n", ...
            solver_str, ...
            idir, ndirection, g, ...
            iseq, nsequence, seq2str(seq), ...
            icmpt, ncompartment, compartments(icmpt));

        % We only keep two points in tlist to output solution at all time
        % steps
        tlist = [0 seq.echotime];

        % Create ODE function and Jacobian from matrices
        [ode_function, Jacobian] = hadc_functions(K{icmpt}, surfint, seq);

        % Set parameters for ODE solver
        options = odeset("Mass", M{icmpt}, "AbsTol", abstol, "RelTol", reltol, ...
            "Vectorized", "on", "Stats", "off", "Jacobian", Jacobian);

        % Solve ODE, keep all time steps (for integral)
        [tvec, y] = solve_ode(ode_function, tlist, rho{icmpt}, options);
        tvec = tvec';
        y = y';

        % Integral over compartment boundary
        hvec = surfint' * y / volumes(icmpt);

        % HADC (free diffusivity minus correction)
        a = trapz(tvec, seq.integral(tvec) .* hvec) / seq.bvalue_no_q;
        adc(iall) = D0 - a;
    end

    % Computational time
    itertimes(iall) = toc(itertime);
end

% Copy results to opposite directions
for idir = dir_inds
    if ~isempty(opposite{idir})
        fprintf("Copying result from direction %d to opposite direction (%d)\n", idir, opposite{idir});
        adc(:, :, opposite{idir}) = adc(:, :, idir);
    end
end

% Compute total HADC (weighted sum over compartments)
weights = initial_density .* volumes; % compartment weights
weights = weights / sum(weights); % normalize
adc_allcmpts(:) = sum(weights' .* adc, 1); % weighted sum

% Create output structure
results.adc = adc;
results.adc_allcmpts = adc_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime);

% Display function evaluation time
toc(starttime);
