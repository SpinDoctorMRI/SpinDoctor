function results = load_hadc(femesh, setup, savepath)
%LOAD_HADC Load the results saved by SOLVE_HADC.
%   The results of each iteration are loaded from
%   "<SAVEPATH>/<GEOMETRYINFO>/<DIFFUSIONINFO>/hadc/<SOLVEROPTIONS>/<SEQUENCEINFO>.MAT".
%
%   setup: struct
%   savepath: string
%   
%   results: struct with fields
%       adc: double(ncompartment, nsequence, ndirection)
%       adc_allcmpts: double(nsequence, ndirection)
%       itertimes: double(ncompartment, nsequence, ndirection)
%       totaltime: double(ncompartment, nsequence, ndirection)


% Measure function evaluation time
starttime = tic;

% Extract experiment parameters
sequences = setup.gradient.sequences;
directions = setup.gradient.directions;
reltol = setup.hadc.reltol;
abstol = setup.hadc.abstol;
solve_ode = setup.hadc.ode_solver;
solver_str = func2str(solve_ode);

initial_density = setup.pde.initial_density;

% Sizes
ncompartment = length(setup.pde.initial_density);
nsequence = length(sequences);
ndirection = size(directions, 2);

% Volumes
volumes = femesh.volumes;

% Folder for loading
savepath = sprintf( ...
    "%s/%s_abstol%g_reltol%g", ...
    savepath, solver_str, abstol, reltol ...
);

% Initialize output arguments
adc = zeros(ncompartment, nsequence, ndirection);
adc_allcmpts = zeros(nsequence, ndirection);
itertimes = zeros(ncompartment, nsequence, ndirection);

inds = [ncompartment ndirection];
for iseq = 1:nsequence
    seq = sequences{iseq};
    % Load results
    filename = sprintf("%s/%s.mat", savepath, seq.string(true));
    fprintf("Load hadc %d/%d.\n", iseq, nsequence);
    mfile = load(filename);
    for iall = 1:prod(inds)
        [icmpt, idir] = ind2sub(inds, iall);
        % Extract parameters for iteration
        ug = directions(:, idir);
        gradient_field = sprintf("cmpt%d_", icmpt) + gradient_fieldstring(ug);
        data = mfile.(gradient_field);

        adc(icmpt, iseq, idir) = data.adc;
        itertimes(icmpt, iseq, idir) = data.itertimes;
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
results.totaltime = sum(itertimes, "all");

% Display function evaluation time
toc(starttime);
