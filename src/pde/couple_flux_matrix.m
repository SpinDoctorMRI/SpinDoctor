function Q = couple_flux_matrix(femesh, pde, Q_blocks, symmetrical)
%COUPLE_FLUX_MATRIX Generate coupling between compartments in flux matrix.
%
%   femesh: struct with fields
%       ncompartment: [1 x 1]
%       point_map: {1 x 4}
%   pde: struct with fields
%       initial_density: [1 x ncompartment]
%       permeability: [1 x nboundary]
%   Q_blocks: {ncompartment x nboundary}[npoint x npoint]
%
%   Q: [ndof_total x ndof_total]


% Mesh quantities
point_map = femesh.point_map;
facets = femesh.facets;

% Sizes
nboundary = femesh.nboundary;
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);
npoint = sum(npoint_cmpts);
cmpt_inds = cumsum([0 npoint_cmpts]);
get_inds = @(icmpt) cmpt_inds(icmpt)+1:cmpt_inds(icmpt+1);

% Create global flux matrix
Q = sparse(npoint, npoint);

% Add boundary contributions to flux matrix
for iboundary = 1:nboundary
    cmpts_touch = find(~cellfun(@isempty, facets(:, iboundary)));
    ntouch = length(cmpts_touch);
    if ntouch == 1
        % Only one compartment touches boundary. It is thus an outer boundary,
        % and may possibly have a boundary relaxation coefficient
        cmpt = cmpts_touch;
        
        % Boundary relaxitity
        k = pde.permeability(iboundary);
        
        % Global indices of the boundary
        inds = get_inds(cmpt);
        
        % Add boundary contribution to global flux matrix for compartment
        Q(inds, inds) = Q(inds, inds) + k * Q_blocks{cmpt, iboundary};
        
    elseif ntouch == 2
        % Two compartments touch boundary, and should be coupled
        cmpt1 = cmpts_touch(1);
        cmpt2 = cmpts_touch(2);
        
        % Initialize the four blocks of the boundary
        Q11 = Q_blocks{cmpt1, iboundary};
        Q12 = sparse(npoint_cmpts(cmpt1), npoint_cmpts(cmpt2));
        % Q21 = sparse(npoint_cmpts(cmpt2), npoint_cmpts(cmpt1));
        Q22 = Q_blocks{cmpt2, iboundary};

        % Identify pairs of points using global numbering
        inds1 = unique(facets{cmpt1, iboundary});
        inds2 = unique(facets{cmpt2, iboundary});
        
        % Check if local indices are already sorted with a one-to-one
        % correspondence
        if all(point_map{cmpt1}(inds1) == point_map{cmpt2}(inds2))
            indinds1 = 1:length(inds1);
            indinds2 = 1:length(inds2);
        else
            % Find coupled nodes
            [indinds1, indinds2] = find(point_map{cmpt1}(inds1) == point_map{cmpt2}(inds2)');
        end
        
        % Create coupling blocks
        Q12(:, inds2(indinds2)) = Q11(:, inds1(indinds1));
        %Q21(:, inds1(indinds1)) = Q22(:, inds2(indinds2));
        Q21 = Q12';
        
        % Check whether to use same permeability coefficients on both sides
        if symmetrical
            % Use the same permeability coefficient on each side of the boundary
            c12 = 1;
            c21 = 1;
        else
            % Weigh permeability coefficients on each side of the boundary with
            % initial spin density equilibrium
            rho1 = pde.initial_density(cmpt1);
            rho2 = pde.initial_density(cmpt2);
            c21 = 2 * rho2 / (rho1 + rho2);
            c12 = 2 * rho1 / (rho1 + rho2);
        end
        
        % Adjust permeability coefficients
        k1 = c21 * pde.permeability(iboundary);
        k2 = c12 * pde.permeability(iboundary);
        
        % Global indices of the boundary in each of the compartments
        inds1 = get_inds(cmpt1);
        inds2 = get_inds(cmpt2);
        
        % Add interface contribution to compartment 1 in global flux matrix
        Q(inds1, inds1) = Q(inds1, inds1) + k1 * Q11;
        Q(inds1, inds2) = Q(inds1, inds2) - k2 * Q12;
        
        % Add interface contribution to compartment 2 in global flux matrix
        Q(inds2, inds1) = Q(inds2, inds1) - k1 * Q21;
        Q(inds2, inds2) = Q(inds2, inds2) + k2 * Q22;
        
        % fprintf("Interface boundary: %d, connecting compartments %d and %d\n", iboundary, cmpt1, cmpt2);
    elseif ntouch > 2
        error("Each interface touch only 1 or 2 compartments");
    end
end
