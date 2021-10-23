function [diffusion_tensor, diffusion_tensor_all, A] = compute_mf_diffusion_tensor(femesh, lap_eig, mf_jn)
%COMPUTE_MF_DIFFUSION_TENSOR Compute the effective diffusion tensor.
%
%   femesh: struct
%   lap_eig: struct
%   mf_jn: cell
%
%   diffusion_tensor: double(3, 3, nsequence, [ncompartment])


% Sizes
nsequence = size(mf_jn{1}, 1);
ncompartment = length(femesh.points);

% Assemble mass matrices
M_cmpts = cell(1, ncompartment);
volumes = zeros(1, ncompartment);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [volumes(icmpt), fevolumes] = get_volume_mesh(points, elements);
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', fevolumes');
end

if length(lap_eig) == 1    % One compartment or some compartments are connected by permeable interfaces
    % Initialize output argument
    diffusion_tensor = zeros(3, 3, nsequence);

    % Assemble mass matrix
    M = blkdiag(M_cmpts{:});
    points = [femesh.points{:}];
    eigfuncs = lap_eig.funcs;
    A = points * M * eigfuncs;

    % Compute diffusion tensor
    for iseq = 1:nsequence
        jn = mf_jn{1}(iseq, :);
        diffusion_tensor(:, :, iseq) = A .* jn * A' / sum(volumes);
    end
    diffusion_tensor_all = diffusion_tensor;
else    % All compartments are uncorrelated
    % Initialize output arguments
    diffusion_tensor = zeros(3, 3, nsequence, ncompartment);
    A = cell(1, ncompartment);

    parfor icmpt = 1:ncompartment
        % Eigenvalues and moments
        M = M_cmpts{icmpt};
        points = femesh.points{icmpt};
        eigfuncs = lap_eig(icmpt).funcs;

        a = points * M * eigfuncs;

        % Compute diffusion tensor
        for iseq = 1:nsequence
            jn = mf_jn{icmpt}(iseq, :);
            diffusion_tensor(:, :, iseq, icmpt) = a .* jn * a' / volumes(icmpt);
        end
        A{icmpt} = a;
    end
    
    volume_fraction = reshape(volumes / sum(volumes), 1, 1, 1, []);
    diffusion_tensor_all = sum(diffusion_tensor.*volume_fraction, 4);
end
