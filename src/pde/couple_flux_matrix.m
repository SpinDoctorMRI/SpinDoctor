function couple_Q = couple_flux_matrix(femesh, boundary_markers, Q)
%COUPLE_FLUX_MATRIX Generate coupling between compartments in flux matrix.
%
%   femesh: struct with fields
%   	ncompartment: int
%    	point_map: cell(1, 4)
%       boundary_markers: cell(ncompartment, nboundary)
% 	Q: cell(1, ncompartment)
%
%   couple_Q: double(ndof_total, ndof_total)


% Mesh quantities
point_map = femesh.point_map;
facets = femesh.facets;

% Sizes
nboundary = femesh.nboundary;
npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);
inds_compartments = cumsum([0 npoint_cmpts]);
get_inds = @(icmpt) inds_compartments(icmpt)+1:inds_compartments(icmpt+1);

% Assemble sparse block diagonal flux matrix
couple_Q = blkdiag(Q{:});

% Couple flux matrix at boundaries
for iboundary = 1:nboundary
    cmpts_touch = find(boundary_markers(:, iboundary));
    ntouch = length(cmpts_touch);
    if ntouch == 2
        % Two compartments touch, and should be coupled
        cmpt1 = cmpts_touch(1);
        cmpt2 = cmpts_touch(2);
        Q1 = Q{cmpt1};
        Q2 = Q{cmpt2};
        Q12 = sparse(size(Q1, 1), size(Q2, 2));
        Q21 = sparse(size(Q2, 1), size(Q1, 2));
        
        % Identify pairs of points using global numbering
        inds1 = unique(facets{cmpt1, iboundary});
        inds2 = unique(facets{cmpt2, iboundary});
        
        if all(point_map{cmpt1}(inds1) == point_map{cmpt2}(inds2))
            indinds1 = 1:length(inds1);
            indinds2 = 1:length(inds2);
        else
            [indinds1, indinds2] = find(point_map{cmpt1}(inds1) == point_map{cmpt2}(inds2)');
        end

        Q12(:, inds2(indinds2)) = -Q1(:, inds1(indinds1));
%         Q21(:, inds1(indinds1)) = -Q2(:, inds2(indinds2));
        
%         % Couple matrices
%         for i = 1:length(inds1)
%             % Loacal indices of coupled points
%             i1 = inds1(i);
%             i2 = inds2(pairs(i, :));
% 
%             % Create flux coupling between points
%             Q12(:, i2) = -Q1(:, i1);
%             % Q21(:, i1) = -Q2(:, i2);
%         end
        
        couple_Q(get_inds(cmpt1), get_inds(cmpt2)) = Q12;
        couple_Q(get_inds(cmpt2), get_inds(cmpt1)) = Q12'; % Q21;
        % fprintf("Interface boundary: %d, connecting compartments %d and %d\n", iboundary, cmpt1, cmpt2);
    elseif ntouch > 2
        error("Each interface touch only 1 or 2 compartments");
    end
end
