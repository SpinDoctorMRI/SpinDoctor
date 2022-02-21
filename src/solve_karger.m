function results = solve_karger(femesh, setup, effective_difftensors)
%SOLVE_KARGER Solve the finite pulse Karger model (FPK).
%
%   femesh: struct
%   setup: struct
%
%   results: struct with fields
%       signal: [ncompartment x namplitude x nsequence x ndirection]
%       signal_allcmpts: [namplitude x nsequence x ndirection]
%       difftensors: [3 x 3 x ncompartment]
%       itertimes: [namplitude x nsequence x ndirection]
%       totaltime: [1 x 1]


% Measure function evaluation time
starttime = tic;

% Extract domain parameters
relaxation = setup.pde.relaxation;
initial_density = setup.pde.initial_density;
permeability = setup.pde.permeability;

% Extract gradient sequences
qvalues = setup.gradient.qvalues;
bvalues = setup.gradient.bvalues;
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;

% Extract experiment parameters
ndirection_karger = setup.karger.ndirection;
reltol = setup.karger.reltol;
abstol = setup.karger.abstol;
solve_ode = setup.karger.ode_solver;
solver_str = func2str(solve_ode);

% Sizes
ncompartment = femesh.ncompartment;
nboundary = femesh.nboundary;
namplitude = setup.namplitude;
nsequence = setup.nsequence;
ndirection = setup.ndirection;

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
A = tau_inv - diag((tau_inv * volumes') ./ volumes');

% Check if user has provided effective diffusion tensors
if nargin == nargin(@solve_karger)
    D = effective_difftensors;
else
    % Compute apparent diffusion coefficients with HADC model
    disp("Computing diffusion tensors using HADC model");
    setup_hadc = setup;
    setup_hadc.hadc.ode_solver = @ode15s;
    setup_hadc.hadc.reltol = 1e-4;
    setup_hadc.hadc.abstol = 1e-6;
    setup_hadc.gradient.sequences = setup_hadc.gradient.sequences(end);
    setup_hadc.gradient.directions = unitsphere(ndirection_karger);
    hadc = solve_hadc(femesh, setup_hadc);

    % Fit effective diffusion tensors
    g = setup_hadc.gradient.directions;
    adc = permute(hadc.adc(:, 1, :), [1 3 2]);
    D = fit_tensor(g, adc);
end

% Relaxation tensor
R = diag(1 ./ relaxation);

% Initial signal
S0 = (volumes .* initial_density).';

% Prepare output structures
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
itertimes = zeros(namplitude, nsequence, ndirection);

% Cartesian indices (for parallel looping with linear indices)
allinds = [namplitude nsequence ndirection];

% Iterate over gradient amplitudes, sequences and directions.
for iall = 1:prod(allinds)
    
    % Measure iteration time
    itertime = tic;

    % Extract Cartesian indices
    [iamp, iseq, idir] = ind2sub(allinds, iall);
    
    % Extract parameters for iteration
    q = qvalues(iamp, iseq);
    b = bvalues(iamp, iseq);
    seq = sequences{iseq};
    g = directions(:, idir);
    
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
        ADC_diag(icmpt, icmpt) = g' * D(:, :, icmpt) * g;
    end
    
    odejac = @(t, y) A - seq.integral(t)^2 * q^2 * ADC_diag - R;
    odefunc = @(t, y) odejac(t, y) * y;
    
    % Set parameters for ODE solver
    options = odeset( ...
        "AbsTol", abstol, ...
        "RelTol", reltol, ...
        "Vectorized", "on", ...
        "Jacobian", odejac, ...
        "Stats", "off");
    
    [~, S] = solve_ode(odefunc, time_list_interval, S0, options);

    signal(:, iall) = S(end, :);
    itertimes(iall) = toc(itertime);
end

% Total parameters
signal_allcmpts(:) = sum(signal, 1);

% Create output structure
results.signal = signal;
results.signal_allcmpts = signal_allcmpts;
results.difftensors = D;
results.itertimes = itertimes;
results.totaltime = toc(starttime);

% Display function evaluation time
toc(starttime);
