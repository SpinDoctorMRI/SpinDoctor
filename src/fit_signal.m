function results = fit_signal(signal, signal_allcmpts, bvalues)
%FIT_SIGNAL Fit the ADC from the dMRI signal
%
%   signal: double(ncompartment, namplitude, nsequence, ndirection)
%   signal_allcmpts: complex double(namplitude, nsequence, ndirection)
%     bvalues: double(namplitude, nsequence)
%
%   results: struct with fields
%       adc: double(ncompartment, nsequence, ndirection)
%          adc_allcmpts_ double(nsequence, ndirection)
%        S0_allcmpts: double(nsequence, ndirection)


% Sizes
[ncompartment, namplitude, nsequence, ndirection] = size(signal);

% Initialize output arguments
adc = zeros(ncompartment, nsequence, ndirection);
adc_allcmpts = zeros(nsequence, ndirection);
S0 = zeros(ncompartment, nsequence, ndirection);
S0_allcmpts = zeros(nsequence, ndirection);

% Only fit ADC if there are two or more b-values
if namplitude == 1
    error("Cannot fit ADC from one b-value only.");
end

for idir = 1:ndirection
    for iseq = 1:nsequence
        b = bvalues(:, iseq)';
        bmin = b(1);
        bmax = b(end);
        for icmpt = 1:ncompartment
            data = real(signal(icmpt, :, iseq, idir));
            [adc_fit, ~, S01d] = process_signal_poly(data, b, bmin, bmax);
            adc(icmpt, iseq, idir) = adc_fit;
            S0(icmpt, iseq, idir) = S01d;
        end
        data = real(signal_allcmpts(:, iseq, idir))';
        [adc_fit, ~, S01d] = process_signal_poly(data, b, bmin, bmax);
        adc_allcmpts(iseq, idir) = adc_fit;
        S0_allcmpts(iseq, idir) = S01d;
    end
end

% Create output structure
results.adc = adc;
results.adc_allcmpts = adc_allcmpts;
results.S0 = S0;
results.S0_allcmpts = S0_allcmpts;
