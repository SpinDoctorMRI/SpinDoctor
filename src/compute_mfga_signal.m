function results = compute_mfga_signal(setup, initial_signal, dtensor)
%COMPUTE_MFGA_SIGNAL Compute the Matrix Formalism Gaussian approximation signal.
%
%   setup: struct with fields
%   initial_density
%   dtensor
%
%   results: struct with fields
%       signal_allcmpts
%       adc_allcmpts


% Extract experiment parameters
bvalues = setup.gradient.bvalues;
directions = setup.gradient.directions;

% Sizes
namplitude = size(bvalues, 1);
nsequence = size(bvalues, 2);
ndirection = size(directions, 2);

% Initialize output arguments
signal = zeros(namplitude, nsequence, ndirection);
adc = zeros(nsequence, ndirection);

% Compute signal
for idir = 1:ndirection
    g = directions(:, idir);
    for iseq = 1:nsequence
        b = setup.gradient.bvalues(:, iseq);
        D = g' * dtensor(:, :, iseq) * g;
        adc(iseq, idir) = D;
        signal(:, iseq, idir) = initial_signal * exp(-D * b);
    end
end

% Create output structure
results.signal_allcmpts = signal;
results.adc_allcmpts = adc;
