function free = compute_free_diffusion(bvalues, diffusivity, volumes, initial_density)
%COMPUTE_FREE_DIFFUSION Compute the free signal and diffusion coefficients.
%
%   bvalues: double(namplitude, nsequence)
% 	diffusivity: double(3, 3, ncompartment)
% 	volumes: double(1, ncompartment)
%  	initial_density: double(1, ncompartment)
%
%   free: struct with fields
%       signal: double(ncompartment, namplitude, nsequence)
%       signal_allcmpts: double(namplitude, nsequence)
%       adc: double(ncompartment, 1)
%       adc_allcmpts: double


% Take direction average diffusivity from tensor (trace)
diffusivity = (diffusivity(1, 1, :) + diffusivity(2, 2, :) + diffusivity(3, 3, :)) / 3;
diffusivity = shiftdim(diffusivity, 1);

% Initial signal
S0 = (initial_density .* volumes).';

% Free signal
signal =  S0 .* exp(-diffusivity' .* shiftdim(bvalues, -1));
signal_allcmpts = shiftdim(sum(signal, 1), 1);

% ADC is simply the diffusion coefficient
adc = diffusivity';

% Total apparent diffusion coefficient weighted by volumes
adc_allcmpts = volumes / sum(volumes) * adc;

% Create output structure
free.signal = signal;
free.signal_allcmpts = signal_allcmpts;
free.adc = adc;
free.adc_allcmpts = adc_allcmpts;
