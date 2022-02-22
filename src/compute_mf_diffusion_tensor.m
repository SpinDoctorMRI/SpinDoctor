function [diffusion_tensor, diffusion_tensor_all, mf_jn, a] = compute_mf_diffusion_tensor(femesh, setup, lap_eig)
%COMPUTE_MF_DIFFUSION_TENSOR Compute the effective diffusion tensor.
%
%   femesh: struct
%   setup: struct
%   lap_eig: struct
%
%   diffusion_tensor: double(3, 3, nsequence, [ncompartment])
%   diffusion_tensor_all: double(3, 3, nsequence)
%   mf_jn: cell(1, ncompartment)
%   A: double(3, neig) or cell(1, ncompartment)


% Compute the JN value that characterizes the contribution of
% a diffusion-encoding sequence to the effective diffusion tensor
mf_jn = compute_mf_jn(lap_eig, setup);

% Get sizes
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
    
    % Compute the matrix A that quantifies the geometrical contribution of
    % a femesh to the effective diffusion tensor
    eigfuncs = lap_eig.funcs;
    a = points * M * eigfuncs;

    % Compute diffusion tensor
    for iseq = 1:nsequence
        jn = mf_jn{1}(iseq, :);
        diffusion_tensor(:, :, iseq) = a .* jn * a' / sum(volumes);
    end
    diffusion_tensor_all = diffusion_tensor;
else    % All compartments are uncorrelated
    % Initialize output arguments
    diffusion_tensor = zeros(3, 3, nsequence, ncompartment);
    a = cell(1, ncompartment);

    % Check if Parallel Computing Toolbox is licensed
    if license('test', 'Distrib_Computing_Toolbox') && isempty(gcp('nocreate'))
        parpool('local', [1, 2048]);
    end

    parfor icmpt = 1:ncompartment
        % Eigenvalues and moments
        M = M_cmpts{icmpt};
        points = femesh.points{icmpt};

        % Compute the matrix 'a' that quantifies the geometrical contribution of
        % a femesh to the effective diffusion tensor
        eigfuncs = lap_eig(icmpt).funcs;
        ai = points * M * eigfuncs;

        % Compute diffusion tensor
        for iseq = 1:nsequence
            jn = mf_jn{icmpt}(iseq, :);
            diffusion_tensor(:, :, iseq, icmpt) = ai .* jn * ai' / volumes(icmpt);
        end
        a{icmpt} = ai;
    end

    volume_fraction = reshape(volumes / sum(volumes), 1, 1, 1, []);
    diffusion_tensor_all = sum(diffusion_tensor.*volume_fraction, 4);
end
