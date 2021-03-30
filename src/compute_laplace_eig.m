function lap_eig = compute_laplace_eig(femesh, pde, eiglim, neig_max)
%COMPUTE_LAPLACE_EIG Compute Laplace eigenvalues, functions and product moments.
%
%   femesh: struct
%   pde: struct
%   eiglim: double
%     neig_max: int
%
%   lap_eig: struct with fields
%       values: cell(1,ndomain)
%       funcs: cell(1,ndomain)
%       moments: cell(1,ndomain)
%       totaltime: double


% Measure computational time of eigendecomposition
starttime = tic;

% Solver parameters
params = {};
% params = [params {"Tolerance" 1e-12}];
params = [params {"Display" true}];

% Check if user has provided a requested number of eigenvalues
if nargin < nargin(@compute_laplace_eig)
    % Compute all eigenvalues
    neig_max = Inf;
end

% Extract domain parameters
diffusivity = pde.diffusivity;
permeability = pde.permeability;
relaxation = pde.relaxation;
boundary_markers = pde.boundary_markers;

% Sizes
ncompartment = femesh.ncompartment;

% Assemble finite element matrices
disp("Setting up FEM matrices");
M_cmpts = cell(1, ncompartment);
K_cmpts = cell(1, ncompartment);
Jx_cmpts = cell(1, 3);
MT2_cmpts = cell(1, ncompartment);
for idim = 1:3
    Jx_cmpts{idim} = cell(1, ncompartment);
end
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    facets = femesh.facets(icmpt, :);
    elements = femesh.elements{icmpt};
    [~, volumes] = get_volume_mesh(points, elements);

    % Assemble flux, stiffness and mass matrices in compartment
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
    K_cmpts{icmpt} = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));

    % Assemble moment matrices (coordinate weighted mass matrices)
    for idim = 1:3
        Jx_cmpts{idim}{icmpt} = mass_matrixP1_3D(elements', volumes', points(idim, :)');
    end

    % Add T2-relaxation term if it is finite and nonnegative. This term is added
    % directly to the stiffness matrix, as it has the effect of adding a
    % decay operator to the diffusion operator
    if isfinite(relaxation(icmpt)) && relaxation(icmpt) > 0
        fprintf("Adding T2-relaxation matrix in compartment %d of %d: %g\n", icmpt, ncompartment, relaxation(icmpt));
        % K_cmpts{icmpt} = K_cmpts{icmpt} + 1 / relaxation(icmpt) * M_cmpts{icmpt};
        MT2_cmpts{icmpt} = 1 / relaxation(icmpt) * M_cmpts{icmpt};
    else
        n = size(M_cmpts{icmpt}, 1);
        MT2_cmpts{icmpt} = sparse(n, n);
    end
end

% Create global mass, stiffness, flux and moment matrices (sparse)
disp("Coupling FEM matrices");
M = blkdiag(M_cmpts{:});
K = blkdiag(K_cmpts{:});
Jx = cellfun(@(J) blkdiag(J{:}), Jx_cmpts, "UniformOutput", false);
MT2 = blkdiag(MT2_cmpts{:});
Q_blocks = assemble_flux_matrix(femesh.points, femesh.facets);
Q = couple_flux_matrix(femesh, pde, Q_blocks, false);

fprintf("Eigendecomposition of FE matrices: size %d x %d\n", size(M));

% % Solve explicit eigenvalue problem, computing all eigenvalues after
% % inverting the mass matrix
% tic
% [funcs, values] = eig(full(M \ (K + Q)));
% toc

% % Solve generalized eigenvalue problem, computing all eigenvalues
% tic
% [funcs, values] = eig(full(K + Q), full(M));
% toc

% Compute at most all eigenvalues in the given domain
neig_max = min(neig_max, size(M, 1));
% ssdim = 4 * neig_max;
% if ssdim < size(M, 1)
%     params = [params {"SubspaceDimension" ssdim}];
% end

% Solve generalized eigenproblem, computing the smallest eigenvalues only.
% If 2 * neig_max >= nnode, a full decomposition is performed,
% calling the eig function inside eigs
tic
[funcs, values] = eigs(K + Q, M, neig_max, "smallestreal", ...
    "IsSymmetricDefinite", true, params{:}); % "smallestabs"
toc

disp("Done with eigendecomposition");

% Order eigenvalues in increasing order
[values, indices] = sort(diag(values));
funcs = funcs(:, indices);

if any(values < 0)
    i = find(values < 0);
    iformat = join(repmat("%d", 1, length(i)));
    vformat = join(repmat("%g", 1, length(i)));
    warning("Found negative eigenvalues: indices " + iformat + ", values " ...
        + vformat + ". Setting them to zero.", ...
        i, values(i));
    values(i) = 0;
end

% Remove eigenvalues above interval defined by length scale
neig_all = length(values);
inds_keep = values <= eiglim;
values = values(inds_keep);
funcs = funcs(:, inds_keep);
neig = length(values);

% Check that the entire interval was explored
if neig == neig_all && ~isinf(eiglim)
    warning("No eigenvalues were outside the interval. Consider increasing neig_max " ...
        + "if there are more eigenvalues that may not have been found in the interval.");
end

fprintf("Found %d eigenvalues on [%g, %g]\n", neig, 0, eiglim);

% Normalize eigenfunctions with mass weighting
disp("Normalizing eigenfunctions");
funcs = funcs ./ sqrt(dot(funcs, M * funcs));

% Compute first order moments of of eigenfunction products
disp("Computing first order moments of products of eigenfunction pairs");
tic
moments = zeros(neig, neig, 3);
for idim = 1:3
    moments(:, :, idim) = funcs' * Jx{idim} * funcs;
end
disp("Computing T2-weighted Laplace mass matrix");
if nnz(MT2)
    massrelax = funcs' * MT2 * funcs;
else
    massrelax = sparse(neig, neig);
end
toc

% Create output structure
lap_eig.values = values;
lap_eig.funcs = funcs;
lap_eig.moments = moments;
lap_eig.massrelax = massrelax;
lap_eig.totaltime = toc(starttime);

% Display function evaluation time
disp("Done with eigendecomposition.");
toc(starttime);
