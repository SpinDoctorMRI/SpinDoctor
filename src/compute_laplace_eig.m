function lap_eig = compute_laplace_eig(femesh, pde, mf, eiglim)
%COMPUTE_LAPLACE_EIG Compute Laplace eigenvalues, functions and product moments.
%
%   femesh: struct
%   pde: struct
%   mf: struct
%   eiglim: [1 x 1]
%
%   lap_eig: struct with fields
%       values: [neig x 1]
%       funcs: [npoint x neig]
%       moments: [neig x neig x 3]
%       massrelax: [neig x neig]
%       totaltime: [1 x 1]


% Measure computational time of eigendecomposition
starttime = tic;

% Solver parameters
params = {};

% Convergence tolerance of eigs
if isfield(mf, 'tolerance') && isnumeric(mf.tolerance)
    params = [params {"Tolerance" mf.tolerance}];
end

% Maximum number of eigs algorithm iterations
if isfield(mf, 'maxiter') && isnumeric(mf.maxiter)
    params = [params {"MaxIterations" mf.maxiter}];
elseif isfield(mf, 'maxiterations') && isnumeric(mf.maxiterations)
    params = [params {"MaxIterations" mf.maxiterations}];
end

% Default: compute all eigenvalues using eig
neig_max = Inf;
% Check if user has provided a requested number of eigenvalues. Finite neig_max triggers eigs.
if isfield(mf, 'neig_max') && isnumeric(mf.neig_max) && mf.neig_max > 0
    neig_max = mf.neig_max;
end
if isinf(neig_max)
    warning("Compute all eigenvalues using EIG which requires much more memory than EIGS. " ...
                + "Eigenvalues out of the interval defined by length scale will be removed.");
end

% Extract domain parameters
diffusivity = pde.diffusivity;
relaxation = pde.relaxation;
no_relaxation = all(isinf(relaxation));

% Sizes
ncompartment = femesh.ncompartment;

% Assemble finite element matrices
disp("Setting up FEM matrices");
M_cmpts = cell(1, ncompartment);
K_cmpts = cell(1, ncompartment);
R_cmpts = cell(1, ncompartment);
Jx_cmpts = repmat({cell(1, ncompartment)}, 1, 3);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [~, volumes] = get_volume_mesh(points, elements);

    % Assemble mass, stiffness, and T2-relaxation matrices in compartment
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
    K_cmpts{icmpt} = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));
    if no_relaxation
        R_cmpts{icmpt} = 0;
    else
        R_cmpts{icmpt} = 1 / relaxation(icmpt) * M_cmpts{icmpt};
    end

    % Assemble moment matrices (coordinate weighted mass matrices)
    for idim = 1:3
        Jx_cmpts{idim}{icmpt} = mass_matrixP1_3D(elements', volumes', points(idim, :)');
    end
end

if all(pde.permeability==0)    % All compartments are uncorrelated
    % Initialize output results
    values_cmpts = cell(1, ncompartment);
    funcs_cmpts = cell(1, ncompartment);
    moments_cmpts = cell(1, ncompartment);
    massrelax_cmpts = cell(1, ncompartment);
    totaltime_cmpts = cell(1, ncompartment);

    parfor icmpt = 1:ncompartment
        starttime_icmpt = tic;

        % Get mass, stiffness, relaxation, flux, and moment matrices (sparse)
        M = M_cmpts{icmpt};
        K = K_cmpts{icmpt};
        R = R_cmpts{icmpt};
        Jx = cellfun(@(J) J{icmpt}, Jx_cmpts, "UniformOutput", false);

        fprintf("Eigendecomposition of FE matrices of %d-th compartment: size %d x %d\n", icmpt, size(M));

        % Compute at most all eigenvalues in the given domain
        params_icmpts = params;
        neig_max_icmpt = min(neig_max, size(M, 1));

        % Solve generalized eigenproblem, computing the smallest eigenvalues only.
        % If 2 * neig_max >= nnode, a full decomposition is performed,
        % calling the eig function inside eigs
        tic
        [funcs, values] = eigs(K, M, neig_max_icmpt, "smallestreal", ...
            "IsSymmetricDefinite", true, params_icmpts{:}); % "smallestabs"
        toc

        disp("Done with eigendecomposition");

        % Order eigenvalues in increasing order
        [values, indices] = sort(diag(values));
        funcs = funcs(:, indices);

        if any(values < 0)
            i = find(values < 0);
            i_str = sprintf(join(repmat("%d", 1, length(i))), i);
            v_str = sprintf(join(repmat("%g", 1, length(i))), values(i));
            warning("Found negative eigenvalues: indices " + i_str + ", values " + v_str + ". Setting them to zero.");
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
        funcs = funcs ./ sqrt(dot(funcs, M * funcs));

        % Compute first order moments of eigenfunction products
        moments = zeros(neig, neig, 3);
        for idim = 1:3
            moments(:, :, idim) = funcs' * Jx{idim} * funcs;
        end
        % Compute T2-weighted Laplace mass matrix
        if no_relaxation
            massrelax = 0;
        else
            massrelax = funcs' * R * funcs;
        end

        % Create output cell
        values_cmpts{icmpt} = values;
        funcs_cmpts{icmpt} = funcs;
        moments_cmpts{icmpt} = moments;
        massrelax_cmpts{icmpt} = massrelax;
        totaltime_cmpts{icmpt} = toc(starttime_icmpt);
    end
    lap_eig = struct('values', values_cmpts, ...
                     'funcs', funcs_cmpts, ...
                     'moments', moments_cmpts, ...
                     'massrelax', massrelax_cmpts, ...
                     'totaltime', totaltime_cmpts);

else    % One compartment or some compartments are connected by permeable interfaces
    % Create global mass, stiffness, relaxation, flux, and moment matrices (sparse)
    disp("Coupling FEM matrices");
    M = blkdiag(M_cmpts{:});
    K = blkdiag(K_cmpts{:});
    Jx = cellfun(@(J) blkdiag(J{:}), Jx_cmpts, "UniformOutput", false);
    Q_blocks = assemble_flux_matrix(femesh.points, femesh.facets);
    Q = couple_flux_matrix(femesh, pde, Q_blocks, false);
    if ~no_relaxation
        R = blkdiag(R_cmpts{:});
    end

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
    if isfield(mf, 'subspacedimension') && isnumeric(mf.subspacedimension)
        params = [params {"SubspaceDimension" mf.subspacedimension}];
    elseif isfield(mf, 'ssdim') && isnumeric(mf.ssdim)
        params = [params {"SubspaceDimension" mf.ssdim}];
    end
    % Display algorithm iteration details
    params = [params {"Display" true}];
    
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
        i_str = sprintf(join(repmat("%d", 1, length(i))), i);
        v_str = sprintf(join(repmat("%g", 1, length(i))), values(i));
        warning("Found negative eigenvalues: indices " + i_str + ", values " + v_str + ". Setting them to zero.");
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
    moments = zeros(neig, neig, 3);
    for idim = 1:3
        moments(:, :, idim) = funcs' * Jx{idim} * funcs;
    end
    disp("Computing T2-weighted Laplace mass matrix");
    if no_relaxation
        massrelax = 0;
    else
        massrelax = funcs' * R * funcs;
    end

    % Create output structure
    lap_eig.values = values;
    lap_eig.funcs = funcs;
    lap_eig.moments = moments;
    lap_eig.massrelax = massrelax;
    lap_eig.totaltime = toc(starttime);
end

% Display function evaluation time
disp("Done with eigendecomposition.");
toc(starttime);
