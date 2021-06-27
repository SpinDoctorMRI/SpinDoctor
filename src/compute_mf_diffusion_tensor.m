function diffusion_tensor = compute_mf_diffusion_tensor(femesh, lap_eig, mf_jn)
%COMPUTE_MF_DIFFUSION_TENSOR Compute the effective diffusion tensor.
%
%   femesh
%   lap_eig
%   mf_jn
%
%   diffusion_tensor: double(3, 3, nsequence)


% Eigenvalues and moments
eigfuncs = lap_eig.funcs;

% Sizes
nsequence = size(mf_jn, 1);
ncompartment = length(femesh.points);

% Initialize output arguments
diffusion_tensor = zeros(3, 3, nsequence);

% Assemble mass matrix
M_cmpts = cell(1, ncompartment);
volumes = zeros(1, ncompartment);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [volumes(icmpt), fevolumes] = get_volume_mesh(points, elements);
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', fevolumes');
end
M = blkdiag(M_cmpts{:});
points = [femesh.points{:}];
a = points * M * eigfuncs;

% Compute diffusion tensor
for iseq = 1:nsequence
    jn = mf_jn(iseq, :);
    diffusion_tensor(:, :, iseq) = a .* jn * a' / sum(volumes);
end
