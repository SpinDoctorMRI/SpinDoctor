function diffusion_tensor = compute_mf_diffusion_tensor(lap_eig, mf_jn, diffusivity)
%COMPUTE_MF_DIFFUSION_TENSOR Compute the effective diffusion tensor.
%
%   lap_eig
%   mf_jn
%   diffusivity
%
%   diffusion_tensor: double(3, 3, nsequence)


% Eigenvalues and moments
moments = lap_eig.moments;

% Sizes
nsequence = size(mf_jn, 1);

% Initialize output arguments
diffusion_tensor = zeros(3, 3, nsequence);

% Compute diffusion tensor
moment = shiftdim(moments(1, :, :), 1);
for iseq = 1:nsequence
    jn = mf_jn(iseq, :);
    diffusion_tensor(:, :, iseq) = diffusivity * moment' .* jn * moment;
end
