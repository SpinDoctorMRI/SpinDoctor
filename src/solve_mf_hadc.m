function mf_hadc = solve_mf_hadc(femesh, setup, lap_eig)
%SOLVE_MF_HADC Solve HADC model using Matrix Formalism.
%
%   femesh: struct
%   setup: struct
%   lap_eig: struct
%
%   mf_hadc: struct
%       signal: double(ncompartment, namplitude, nsequence, ndirection)
%            compartment-wise mfga signal.
%       signal_allcmpts: double(namplitude, nsequence, ndirection)
%           mfga signal.
%       adc: double(ncompartment, nsequence, ndirection)
%           compartment-wise adc.
%       adc_allcmpts: double(nsequence, ndirection)
%           adc.
%       diffusion_tensor: double(3, 3, nsequence, [ncompartment])
%           compartment-wise diffusion tensor.
%       diffusion_tensor_allcmpts: double(3, 3, nsequence)
%           diffusion tensor.
%       jn: cell(1, ncompartment)
%           jn of adc.
%       a: double(3, neig) or cell(1, ncompartment)
%           moment a of adc.

% TEMPORARY. Camino file sequences not yet implemented for this solver.
const_ind = cellfun(@(x) ~isa(x,"SequenceCamino"),setup.gradient.sequences,'UniformOutput',true);
if ~all(const_ind,'all')
    warning("Currently %s does not support camino file sequences. \n Solving only for non-camino sequences",mfilename);
    setup.gradient.sequences = setup.gradient.sequences(const_ind);
    setup.nsequence = sum(const_ind);
end

% Compute the Matrix Formalism effective diffusion tensor, jn and moment a
[diffusion_tensor, diffusion_tensor_all, mf_jn, a] = ...
            compute_mf_diffusion_tensor(femesh, setup, lap_eig);

% Extract experiment parameters
bvalues = setup.gradient.bvalues;
directions = setup.gradient.directions;
initial_signal = setup.pde.initial_signal;

% Sizes
ncompartment = setup.ncompartment;
namplitude = setup.namplitude;
nsequence = setup.nsequence;
ndirection = setup.ndirection;

% Initialize output arguments
signal = zeros(ncompartment, namplitude, nsequence, ndirection);
signal_allcmpts = zeros(namplitude, nsequence, ndirection);
adc = zeros(ncompartment, nsequence, ndirection);
adc_allcmpts = zeros(nsequence, ndirection);

% Compute signal_allcmpts
for idir = 1:ndirection
    g = directions(:, idir);
    for iseq = 1:nsequence
        b = bvalues(:, iseq);
        D = g' * diffusion_tensor_all(:, :, iseq) * g;
        adc_allcmpts(iseq, idir) = D;
        signal_allcmpts(:, iseq, idir) = initial_signal * exp(-D * b);
    end
end

% Compute compartmentwise signal
if length(size(diffusion_tensor)) == 4
    for idir = 1:ndirection
        g = directions(:, idir);
        for iseq = 1:nsequence
            b = bvalues(:, iseq);
            for icmpt = 1:ncompartment
                D = g' * diffusion_tensor(:, :, iseq, icmpt) * g;
                adc(icmpt, iseq, idir) = D;
                signal(icmpt, :, iseq, idir) = initial_signal * exp(-D * b);
            end
        end
    end
else
    adc = adc_allcmpts;
    signal = signal_allcmpts;
end

% Create output structure
mf_hadc.signal = signal;
mf_hadc.signal_allcmpts = signal_allcmpts;
mf_hadc.adc = adc;
mf_hadc.adc_allcmpts = adc_allcmpts;
mf_hadc.diffusion_tensor = diffusion_tensor;
mf_hadc.diffusion_tensor_allcmpts = diffusion_tensor_all;
mf_hadc.jn = mf_jn;
mf_hadc.a = a;
end
