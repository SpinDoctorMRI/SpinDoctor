function [adc, adc_allcmpts] = compute_adc_sta(femesh, setup)
%COMPUTE_ADC_STA Compute the ADC in the short diffusion time regime.
%
%   femesh: struct
%   setup: struct
%
%   adc: double(ncompartment, nsequence, ndirection)
%     adc_allcmpts: double(nsequence, ndirection)


% Extract domain parameters
diffusivity = setup.pde.diffusivity;
initial_density = setup.pde.initial_density;

% Exeriment paramters
sequences = setup.gradient.sequences;

% Gradient directions
dir_points = setup.gradient.directions.points;

% Sizes
nsequence = length(sequences);
ncompartment = femesh.ncompartment;
ndirection = setup.gradient.ndirection;

% Prepare output structures
adc = zeros(ncompartment, nsequence, ndirection);
adc_allcmpts = zeros(nsequence, ndirection);

% Compute ADC
volumes = zeros(1, ncompartment);
for icmpt = 1:ncompartment
    % FE Mesh
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    facets = [femesh.facets{icmpt, :}];

    % Get total volume and surface areas with normals for each facet
    volumes(icmpt) = get_volume_mesh(points, elements);
    [~, areas, ~, normals] = get_surface_mesh(points, facets);

    for idir = 1:ndirection
        
        g = dir_points(:, idir);
        
        D0 = g' * diffusivity(:, :, icmpt) * g;

        % Project surface areas onto plane orthogonal to gradient direction
        SAu = (g' * normals).^2 * areas';

        for iseq = 1:nsequence
            seq = sequences{iseq};
            delta = seq.delta;
            Delta = seq.Delta;
            if isa(seq, "PGSE")
                D = D0 - D0^(3/2) * SAu / volumes(icmpt) * 16 / 35 / sqrt(pi) ...
                    / delta^2 / (3 * Delta - delta) * ( ...
                      (Delta - delta)^(7/2) ...
                    + (Delta + delta)^(7/2) ...
                    - 2 * (delta^(7/2) + Delta^(7/2)));
            else
                diffusion_time = seq.diffusion_time;
                D = D0 * (1 - 4 / 3 * sqrt(D0 / pi) * SAu / volumes(icmpt) * sqrt(diffusion_time));
            end
            adc(icmpt, iseq, idir) = D;
        end
    end
end

% Compute initial signal weigthed global ADC from compartments
weights = (initial_density .* volumes)';
weights = weights / sum(weights);
adc_allcmpts(:, :) = sum(weights .* adc, 1);
