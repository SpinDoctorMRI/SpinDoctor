function lap_eig = compute_laplace_eig(femesh, pde, mf, savepath)
%COMPUTE_LAPLACE_EIG Compute Laplace eigenvalues and eigenfunctions.
%
%   femesh: struct
%   pde: struct
%   mf: struct
%   savepath: path string
%
%   lap_eig: struct with fields
%       values: [neig x 1]
%       funcs: [npoint x neig]
%       length_scales: [neig x 1]
%       totaltime: [1 x 1]


% Measure computational time of eigendecomposition
starttime = tic;

% Check if a save path has been provided (this toggers saving)
do_save = nargin == nargin(@compute_laplace_eig);

% Extract mf parameters
rerun = mf.rerun_eigen;
if mf.surf_relaxation
    % surface relaxation is on, use lap_eig with zero Neumann condition
    pde.permeability = 0;
end

% Solver parameters
params = {};

% set eigs' sigma
if isfield(mf, 'eigs')
    params = [params {mf.eigs.sigma}];
else
    params = [params {-1e-8}];
end

% Display algorithm iteration details
params = [params {"Display" true}];

if isfield(mf, 'eigs')
    % Convergence tolerance of eigs
    if isfield(mf.eigs, 'tolerance') && isnumeric(mf.eigs.tolerance)
        params = [params {"Tolerance" mf.eigs.tolerance}];
    end
    % Maximum number of eigs algorithm iterations
    if isfield(mf.eigs, 'maxiter') && isnumeric(mf.eigs.maxiter)
        params = [params {"MaxIterations" mf.eigs.maxiter}];
    elseif isfield(mf.eigs, 'maxiterations') && isnumeric(mf.eigs.maxiterations)
        params = [params {"MaxIterations" mf.eigs.maxiterations}];
    end
end

% Extract domain parameters
diffusivity = pde.diffusivity;
mean_diffusivity = pde.mean_diffusivity;

% Compute the upper bound of the eigenvalue interval
eiglim = length2eig(mf.length_scale, mean_diffusivity);

% Sizes
ncompartment = femesh.ncompartment;

% Check if saved lap_eig exists
no_result = true;
if do_save && ~rerun
    lap_eig = load_laplace_eig(savepath, mf, mean_diffusivity);
    if ~isempty(lap_eig)
        no_result = false;
    end
end

if no_result
    if isinf(mf.neig_max)
        warning("Computing all eigenvalues using EIG requires much more memory than EIGS." ...
                + " Eigenvalues out of the interval defined by length scale will be removed.");
    end

    % Assemble finite element matrices
    disp("Setting up FEM matrices");
    M_cmpts = cell(1, ncompartment);
    K_cmpts = cell(1, ncompartment);
    for icmpt = 1:ncompartment
        % Finite elements
        points = femesh.points{icmpt};
        elements = femesh.elements{icmpt};
        [~, volumes] = get_volume_mesh(points, elements);

        % Assemble mass and stiffness matrices in compartment
        M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
        K_cmpts{icmpt} = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));
    end

    if all(pde.permeability==0)    % All compartments are uncorrelated
        % Initialize output results
        lap_eig = struct();

        for icmpt = 1:ncompartment
            starttime_icmpt = tic;

            % Get mass and stiffness matrices (sparse)
            M = M_cmpts{icmpt};
            K = K_cmpts{icmpt};

            fprintf("Eigendecomposition of FE matrices of compartment %d/%d: size %d x %d\n", ...
                icmpt, ncompartment, size(M));

            [eigvals, eigfuncs] = eigendecompose(K, M, params, mf, eiglim);

            % Save eigendecomposition
            lap_eig(icmpt).values = eigvals;
            lap_eig(icmpt).funcs = eigfuncs;
            lap_eig(icmpt).totaltime = toc(starttime_icmpt);
        end

    else    % One compartment or some compartments are connected by permeable interfaces
        % Create global mass, stiffness, and flux matrices (sparse)
        disp("Coupling FEM matrices");
        M = blkdiag(M_cmpts{:});
        K = blkdiag(K_cmpts{:});
        Q_blocks = assemble_flux_matrix(femesh.points, femesh.facets);
        Q = couple_flux_matrix(femesh, pde, Q_blocks, false);

        fprintf("Eigendecomposition of FE matrices: size %d x %d\n", size(M));

        [eigvals, eigfuncs] = eigendecompose(K + Q, M, params, mf, eiglim);

        % Create output structure
        lap_eig.values = eigvals;
        lap_eig.funcs = eigfuncs;
        lap_eig.totaltime = toc(starttime);
    end
    
    % add length scales
    lap_eig = add_eig_length(lap_eig, mean_diffusivity);
end

% get actual length scale
if isinf(mf.neig_max)
    ls = mf.length_scale;
else
    ls_min = cellfun(@(x) min(x), {lap_eig(:).length_scales});
    ls = max(ls_min);
end

% save laplace eigendecomposition
if do_save && no_result
    if ~isdir(savepath)
        warning('%s does not exist. Creating it to store eigenfunctions\n',savepath);
        mkdir(savepath)
    end
    % file for saving
    filename = sprintf( ...
        "%s/lap_eig_lengthscale%.4f_neigmax%g.mat", ...
        savepath, ls, mf.neig_max ...
    );

    disp("save " + filename + " -v7.3 lap_eig");
    save(filename, "-v7.3", "lap_eig");
end

% reset lap_eig according to eiglim
if ls < mf.length_scale
    disp("Remove eigenvalues above interval defined by length scale.");
    lap_eig = reset_lapeig(lap_eig, eiglim, Inf);
end

% Display function evaluation time
disp("Done with eigendecomposition.");
toc(starttime);
end

function [values, funcs] = eigendecompose(KQ, M, params, mf, eiglim)
%EIGENDECOMPOSE Solve the generalized eigenvalue problem KQ*V = M*V*D.
%
%   KQ: matrix [npoint x npoint]
%   M: matrix [npoint x npoint]
%   params: cell - parameters of eigs
%   mf: struct
%   eiglim: double - maximal eigenvalues
%
%   values: [neig x 1]
%   funcs: [npoint x neig]


% Set number of eigenvalues to compute
neig_max = min([mf.neig_max, size(M, 1)]);
% Set maximum size of Krylov subspace
if isfield(mf, 'eigs') && isfield(mf.eigs, 'subspacedimension') ...
    && isnumeric(mf.eigs.subspacedimension)

    ssdim = max([neig_max+2, mf.eigs.subspacedimension]);
    if ssdim < size(M, 1)
        params = [params {"SubspaceDimension" ssdim}];
    else
        warning("SubspaceDimension is invalid (> dim of mass matrix).")
    end
elseif isfield(mf, 'eigs') && isfield(mf.eigs, 'ssdim') ...
    && isnumeric(mf.eigs.ssdim)

    ssdim = max([neig_max+2, mf.eigs.ssdim]);
    if ssdim < size(M, 1)
        params = [params {"SubspaceDimension" ssdim}];
    else
        warning("SubspaceDimension is invalid (> dim of mass matrix).")
    end
end

% Solve generalized eigenproblem, computing the smallest eigenvalues only.
% If 2 * neig_max >= nnode, a full decomposition is performed,
% calling the eig function inside eigs
tic
[funcs, values] = eigs(KQ, M, neig_max, params{:}, ...
    "IsSymmetricDefinite", true);
toc

% Sort eigenvalues in increasing order
[values, indices] = sort(diag(values));
funcs = funcs(:, indices);

if any(values < 0)
    i = find(values < 0);
    i_str = sprintf(join(repmat("%d", 1, length(i))), i);
    v_str = sprintf(join(repmat("%g", 1, length(i))), values(i));
    warning("Found negative eigenvalues: indices " + i_str + ", values " + v_str + ". Setting them to zero.");
    values(i) = 0;
end

if isinf(mf.neig_max)
    % Remove eigenvalues above interval defined by length scale
    inds_keep = values <= eiglim;
    values = values(inds_keep);
    funcs = funcs(:, inds_keep);
    neig = length(values);
else
    % keep all eigenvalues if mf.neig_max is not inf
    inds_keep = values <= eiglim;
    neig = length(values);
end

% Check that the entire interval was explored
if all(inds_keep) && ~isinf(eiglim)
    warning("No eigenvalues were outside the interval. Consider increasing neig_max " ...
        + "if there are more eigenvalues that may not have been found in the interval.");
end

fprintf("Found %d eigenvalues on [%g, %g].\n", neig, 0, max(values));

% Normalize eigenfunctions with mass weighting
funcs = funcs ./ sqrt(dot(funcs, M * funcs));
end
