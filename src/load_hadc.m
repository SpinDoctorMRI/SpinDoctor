function results = load_hadc(femesh, setup, savepath)
%LOAD_HADC Load the results saved by SOLVE_HADC.
%   The results of each iteration are loaded from
%   "<SAVEPATH>/HADC_<SOLVER_OPTIONS>/<ITERATION_INFO>.MAT".
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

initial_density = setup.pde.initial_density;

% Sizes
ncompartment = length(setup.pde.initial_density);
nsequence = length(sequences);
ndirection = size(directions, 2);

% Volumes
volumes = zeros(1, ncompartment);
for icmpt = 1:ncompartment
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    volumes(icmpt) = get_volume_mesh(points, elements);
end

% Folder for loading
savepath = sprintf( ...
    "%s/hadc_abstol%g_reltol%g", ...
    savepath, abstol, reltol ...
);
    
% Initialize output arguments
adc = zeros(ncompartment, nsequence, ndirection);
adc_allcmpts = zeros(nsequence, ndirection);
itertimes = zeros(ncompartment, nsequence, ndirection);

% Cartesian indices (for parallel looping with linear indices)
allinds = [ncompartment nsequence ndirection];

% Iterate over gradient amplitudes, sequences and directions. If the Matlab
% PARALLEL COMPUTING TOOLBOX is available, the iterations may be done in
% parallel, otherwise it should work like a normal loop. If that is not the
% case, replace the `parfor` keyword by the normal `for` keyword.
parfor iall = 1:prod(allinds)

    % Extract indices
    [icmpt, iseq, idir] = ind2sub(allinds, iall);

    % Extract parameters for iteration
    seq = sequences{iseq};
    g = directions(:, idir);

    % File name for saving or loading iteration results
    filename = sprintf("%s/%s.mat", savepath, adc_string(icmpt, seq, g));
    
    % Load results
    fprintf("Load %s\n", filename);
    mfile = matfile(filename, "Writable", false);
    adc(iall) = mfile.adc;
    itertimes(iall) = mfile.itertime;
end % iterations

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
