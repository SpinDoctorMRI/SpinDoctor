function lap_eig = reset_lapeig(femesh, pde, lap_eig, eiglim)
%RESET_LAPEIG Reset Laplace eigenvalues, functions and product moments based on a new eiglim.
%
%   femesh: struct
%   pde: struct
%   lap_eig: struct
%   eiglim: [1 x 1]
%
%   lap_eig: struct with fields
%       values: [neig x 1]
%       funcs: [npoint x neig]
%       moments: [neig x neig x 3]
%       massrelax: [neig x neig] or 0
%       totaltime: [1 x 1]
%       length_scale (optional): [neig x 1]


% Measure computational time of eigendecomposition
starttime = tic;

% Extract domain parameters
relaxation = pde.relaxation;
no_relaxation = all(isinf(relaxation));

% Sizes
ncompartment = femesh.ncompartment;

% Assemble finite element matrices
M_cmpts = cell(1, ncompartment);
R_cmpts = cell(1, ncompartment);
Jx_cmpts = repmat({cell(1, ncompartment)}, 1, 3);
for icmpt = 1:ncompartment
    % Finite elements
    points = femesh.points{icmpt};
    elements = femesh.elements{icmpt};
    [~, volumes] = get_volume_mesh(points, elements);

    % Assemble mass, stiffness, and T2-relaxation matrices in compartment
    M_cmpts{icmpt} = mass_matrixP1_3D(elements', volumes');
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
    for icmpt = 1:ncompartment
        % Get mass, stiffness, relaxation, flux, and moment matrices (sparse)
        R = R_cmpts{icmpt};
        Jx = cellfun(@(J) J{icmpt}, Jx_cmpts, "UniformOutput", false);

        % Load eigenvalues and eigenfunctions
        values = lap_eig(icmpt).values;
        funcs = lap_eig(icmpt).funcs;

        % Remove eigenvalues above interval defined by length scale
        neig_all = length(values);
        inds_keep = values <= eiglim;

        if all(inds_keep, 'all')
            warning('Compartment %d: no eigenvalue is removed.', icmpt);
        else
            % Remove moments and massrelax
            lap_eig(icmpt).moments = 0;
            lap_eig(icmpt).massrelax = 0;

            % Remove out-of-range eigenvalues
            values = values(inds_keep);
            funcs = funcs(:, inds_keep);
            neig = length(values);
            if isfield(lap_eig, 'length_scales')
                lap_eig(icmpt).length_scales = lap_eig(icmpt).length_scales(inds_keep);
            end
            fprintf("Compartment %d: remove %d eigenvalues.\n", icmpt, neig_all - neig);
            
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
            
            % Reset lap_eig
            lap_eig(icmpt).values = values;
            lap_eig(icmpt).funcs = funcs;
            lap_eig(icmpt).moments = moments;
            lap_eig(icmpt).massrelax = massrelax;
        end
    end
else    % One compartment or some compartments are connected by permeable interfaces
    % Create global mass, stiffness, relaxation, flux, and moment matrices (sparse)
    Jx = cellfun(@(J) blkdiag(J{:}), Jx_cmpts, "UniformOutput", false);
    if ~no_relaxation
        R = blkdiag(R_cmpts{:});
    end

    % Load eigenvalues and eigenfunctions
    values = lap_eig.values;
    funcs = lap_eig.funcs;

    % Remove eigenvalues above interval defined by length scale
    neig_all = length(values);
    inds_keep = values <= eiglim;

    if all(inds_keep, 'all')
        warning('No eigenvalue is removed.');
    else
        % Remove moments and massrelax
        lap_eig.moments = 0;
        lap_eig.massrelax = 0;

        % Remove out-of-range eigenvalues
        values = values(inds_keep);
        funcs = funcs(:, inds_keep);
        neig = length(values);
        if isfield(lap_eig, 'length_scales')
            lap_eig.length_scales = lap_eig.length_scales(inds_keep);
        end
        fprintf("Remove %d eigenvalues.\n", neig_all - neig);
        
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

        % Reset lap_eig
        lap_eig.values = values;
        lap_eig.funcs = funcs;
        lap_eig.moments = moments;
        lap_eig.massrelax = massrelax;
    end
end

% Display function evaluation time
disp("Done with reset.");
toc(starttime);
