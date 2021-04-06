function results = solve_karger(femesh, setup)
%SOLVE_KARGER Solve the finite pulse Karger model.
%
%   femesh: struct
%   setup: struct
%
%   results: struct with fields
%       signal: [ncompartment x namplitude x nsequence x ndirection]
%       signal_allcmpts: [namplitude x nsequence x ndirection]
%       difftensors: [3 x 3 x ncompartment x nsequence x ndirection]
%       difftensors_allcmpts: [3 x 3 x nsequence x ndirection]
%       itertimes: [namplitude x nsequence x ndirection]
%       totaltime: [1 x 1]


% Measure function evaluation time
starttime = tic;

% Only works for strait cylinders for now
assert(setup.geometry.cell_shape == "cylinder");
assert(~any(setup.geometry.deformation));

% Extract domain parameters
diffusivity = setup.pde.diffusivity;
relaxation = setup.pde.relaxation;
initial_density = setup.pde.initial_density;
permeability = setup.pde.permeability;
boundary_markers = setup.pde.boundary_markers;

% Extract gradient sequences
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;

% Extract HARDI directions
dir_points = setup.gradient.directions.points;
dir_inds = setup.gradient.directions.indices;
opposite = setup.gradient.directions.opposite;

% Extract experiment parameters
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
reltol = setup.karger.reltol;
abstol = setup.karger.abstol;
solve_ode = setup.karger.ode_solver;
solver_str = func2str(solve_ode);

% Sizes
ncompartment = femesh.ncompartment;
nboundary = femesh.nboundary;
namplitude = size(qvalues, 1);
nsequence = length(sequences);
ndirection = setup.gradient.ndirection;
ndirection_unique = length(dir_inds);

% Compute volumes and surface areas
volumes = zeros(1, ncompartment);
surface_areas = sparse(ncompartment, nboundary);
for icmpt = 1:ncompartment
    % FE Mesh
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    facets = femesh.facets(icmpt, :);

    % Volume
    volumes(icmpt) = get_volume_mesh(points, elements);

    % Surface area
    surface_areas(icmpt, :) = cellfun(@(f) get_surface_mesh(points, f), facets);
end

% Associate surface areas to the compartment interfaces they measure
tau_inv = sparse(ncompartment, ncompartment);
for iboundary = 1:nboundary
    inds = find(surface_areas(:, iboundary));
    if length(inds) == 2
        i1 = inds(1);
        i2 = inds(2);
        tmp = permeability(iboundary) * surface_areas(i1, iboundary);
        tau_inv(i1, i2) = tmp;
        tau_inv(i2, i1) = tmp;
    end
end
tau_inv = tau_inv ./ volumes;
A = tau_inv - diag(sum(volumes .* tau_inv, 2) ./ volumes');

% Compute effective diffusion tensors in compartments
Deff_cmpts = zeros(3, 3, ncompartment);
for icmpt = 1:ncompartment
%     % Finite elements
%     points = femesh.points{icmpt};
%     facets = femesh.facets(icmpt, :);
%     elements = femesh.elements{icmpt};
%     [~, fe_volumes] = get_volume_mesh(points, elements);
% 
%     % Assemble flux, stiffness and mass matrices in compartment
%     M = mass_matrixP1_3D(elements', fe_volumes');
%     K = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));
%     % Q = assemble_flux_matrix(points, facets, permeability);

    % Deff_cmpts(3, 3, icmpt) = diffusivity(3, 3, icmpt);
end

% Compute apparent diffusion coefficient with HADC model
setup_hadc = setup;
setup_hadc.hadc.ode_solver = @ode15s;
setup_hadc.hadc.reltol = 1e-4;
setup_hadc.hadc.abstol = 1e-4;
setup_hadc.gradient.sequences = setup_hadc.gradient.sequences(1);
setup_hadc.gradient.directions.points = [
    1 0 0 1 1 0
    0 1 0 1 0 1
    0 0 1 0 1 1
];
setup_hadc.gradient.ndirection = 6;
setup_hadc.gradient.directions.indices = 1:6;
setup_hadc.gradient.directions.opposite = cell(1, 6);
hadc = solve_hadc(femesh, setup_hadc);

% Deduce effective diffusion tensor in each compartment
a = permute(hadc.adc(:, 1, :), [1 3 2]);
Deff_cmpts(1, 1, :) = a(:, 1);
Deff_cmpts(2, 2, :) = a(:, 2);
Deff_cmpts(3, 3, :) = a(:, 3);
Deff_cmpts(1, 2, :) = (a(:, 4) - a(:, 1) - a(:, 2)) / 2;
Deff_cmpts(2, 1, :) = (a(:, 4) - a(:, 1) - a(:, 2)) / 2;
Deff_cmpts(1, 3, :) = (a(:, 4) - a(:, 1) - a(:, 3)) / 2;
Deff_cmpts(3, 1, :) = (a(:, 4) - a(:, 1) - a(:, 3)) / 2;
Deff_cmpts(2, 3, :) = (a(:, 4) - a(:, 2) - a(:, 3)) / 2;
Deff_cmpts(3, 2, :) = (a(:, 4) - a(:, 2) - a(:, 3)) / 2;

% Prepare output structures
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

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
    
    % Extract parameters for iteration
    q = qvalues(iamp, iseq);
    b = bvalues(iamp, iseq);
    seq = sequences{iseq};
    g = dir_points(:, idir);
    
    % Display state of iterations
    fprintf("Computing Karger signal using %s\n" ...
        + "  Direction %d of %d: g = [%.2f; %.2f; %.2f]\n" ...
        + "  Sequence  %d of %d: f = %s\n" ...
        + "  Amplitude %d of %d: q = %g, b = %g\n", ...
        solver_str, ...
        idir, ndirection, g, ...
        iseq, nsequence, seq, ...
        iamp, namplitude, q, b);
    
    TE = seq.echotime;
    time_list_interval = [0 TE/2 TE];
    
    ADC_diag = sparse(ncompartment, ncompartment);
    for icmpt = 1:ncompartment
        ADC_diag(icmpt, icmpt) = g' * Deff_cmpts(:, :, icmpt) * g;
        % ADC_diag(icmpt, icmpt) = g' * diffusivity(:, :, icmpt) * g;
    end
    
    odejac = @(t, y) A - seq.integral(t) * q^2 * ADC_diag;
    odefunc = @(t, y) odejac(t, y) * y;
    
    % Set parameters for ODE solver
    options = odeset( ...
        "AbsTol", abstol, ...
        "RelTol", reltol, ...
        "Vectorized", "on", ...
        "Jacobian", odejac, ...
        "Stats", "off");
    
    % Initial signal
    S0 = (volumes .* initial_density).';
    
    [~, S] = solve_ode(odefunc, time_list_interval, S0, options);

    signal(:, iall) = S(end, :);
end

% Total parameters
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.itertimes = itertimes;
results.totaltime = toc(starttime);

% Display function evaluation time
toc(starttime);
