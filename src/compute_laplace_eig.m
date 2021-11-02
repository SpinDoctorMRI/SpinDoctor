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
%       md5: md5 hash string


% Measure computational time of eigendecomposition
starttime = tic;
disp("Computing or loading the Laplace eigenfunctions");

% Check if a save path has been provided (this toggers saving)
do_save = nargin == nargin(@solve_btpde);

% Provide default value
rerun = mf.rerun_eigen;

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

% Extract mf parameters
% Set neig_max to Inf to trigger the full eigendecomposition.
neig_max = mf.neig_max;
length_scale = mf.length_scale;
if isinf(neig_max)
    warning("Compute all eigenvalues using EIG which requires much more memory than EIGS. " ...
                + "Eigenvalues out of the interval defined by length scale will be removed.");
end

% Extract domain parameters
diffusivity = pde.diffusivity;
mean_diffusivity = pde.mean_diffusivity;

% Compute the upper bound of the eigenvalue interval
eiglim = length2eig(length_scale, mean_diffusivity);

% Sizes
ncompartment = femesh.ncompartment;

no_result = true;
if do_save
    % Folder for saving
    filename = sprintf( ...
        "%s/lap_eig_lengthscale%g_neigmax%g.mat", ...
        savepath, length_scale, neig_max ...
    );

    if ~rerun
        lap_eig = load_laplace_eig(savepath, mf, mean_diffusivity);
        if ~isempty(lap_eig)
            no_result = false;
        end
    end
end

if no_result
    % Assemble finite element matrices
    disp("Setting up FEM matrices");
    M_cmpts = cell(1, ncompartment);
    K_cmpts = cell(1, ncompartment);
    for icmpt = 1:ncompartment
        % Finite elements
        points = femesh.points{icmpt};
        elements = femesh.elements{icmpt};
        [~, volumes] = get_volume_mesh(points, elements);

        % Assemble mass, stiffness, and T2-relaxation matrices in compartment
        M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
        K_cmpts{icmpt} = stiffness_matrixP1_3D(elements', points', diffusivity(:, :, icmpt));
    end

    if all(pde.permeability==0)    % All compartments are uncorrelated
        % Initialize output results
        values_cmpts = cell(1, ncompartment);
        funcs_cmpts = cell(1, ncompartment);
        totaltime_cmpts = cell(1, ncompartment);
        md5_cmpts = cell(1, ncompartment);

        parfor icmpt = 1:ncompartment
            starttime_icmpt = tic;

            % Get mass, stiffness, relaxation, flux, and moment matrices (sparse)
            M = M_cmpts{icmpt};
            K = K_cmpts{icmpt};

            fprintf("Eigendecomposition of FE matrices of %d-th compartment: size %d x %d\n", icmpt, size(M));

            % Compute at most all eigenvalues in the given domain
            params_icmpts = params;
            neig_max_icmpt = min(neig_max, size(M, 1));

            % Solve generalized eigenproblem, computing the smallest eigenvalues only.
            % If 2 * neig_max >= nnode, a full decomposition is performed,
            % calling the eig function inside eigs
            tic
            [funcs, values] = eigs(K, M, neig_max_icmpt, "smallestreal", ...
                "IsSymmetricDefinite", true, params_icmpts{:});
            toc

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

            % Create output cell
            values_cmpts{icmpt} = values;
            funcs_cmpts{icmpt} = funcs;
            totaltime_cmpts{icmpt} = toc(starttime_icmpt);
            md5_cmpts{icmpt} = DataHash(funcs);
        end
        lap_eig = struct('values', values_cmpts, ...
                        'funcs', funcs_cmpts, ...
                        'totaltime', totaltime_cmpts, ...
                        'md5', md5_cmpts);

    else    % One compartment or some compartments are connected by permeable interfaces
        % Create global mass, stiffness, and flux matrices (sparse)
        disp("Coupling FEM matrices");
        M = blkdiag(M_cmpts{:});
        K = blkdiag(K_cmpts{:});
        Q_blocks = assemble_flux_matrix(femesh.points, femesh.facets);
        Q = couple_flux_matrix(femesh, pde, Q_blocks, false);

        fprintf("Eigendecomposition of FE matrices: size %d x %d\n", size(M));

        % Compute at most all eigenvalues in the given domain
        if isfield(mf, 'subspacedimension') && isnumeric(mf.subspacedimension) && ~isinf(neig_max)
            ssdim = max([neig_max+2, mf.subspacedimension]);
            ssdim = min([ssdim, size(M, 1)]);
            params = [params {"SubspaceDimension" ssdim}];
        elseif isfield(mf, 'ssdim') && isnumeric(mf.ssdim) && ~isinf(neig_max)
            ssdim = max([neig_max+2, mf.ssdim]);
            ssdim = min([ssdim, size(M, 1)]);
            params = [params {"SubspaceDimension" ssdim}];
        end
        neig_max = min([neig_max, size(M, 1)]);
        
        % Display algorithm iteration details
        params = [params {"Display" true}];
        
        % Solve generalized eigenproblem, computing the smallest eigenvalues only.
        % If 2 * neig_max >= nnode, a full decomposition is performed,
        % calling the eig function inside eigs
        tic
        [funcs, values] = eigs(K + Q, M, neig_max, "smallestreal", ...
            "IsSymmetricDefinite", true, params{:});
        toc

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

        % Create output structure
        lap_eig.values = values;
        lap_eig.funcs = funcs;
        lap_eig.totaltime = toc(starttime);
        lap_eig.md5 = DataHash(funcs);
    end
end

% add length scales
lap_eig = add_eig_length(lap_eig, mean_diffusivity);

% save laplace eigendecomposition
if do_save && no_result
    disp("save " + filename + " -v7.3 lap_eig");
    save(filename, "-v7.3", "lap_eig");
end

% Display function evaluation time
disp("Done with eigendecomposition.");
toc(starttime);
