function flux_matrix = assemble_flux_matrix(points, facets, permeability)
%ASSEMBLE_FLUX_MATRIX Assemble flux matrix in compartment
%
%   points: double(3, nnode)
%   facets: cell(1, nboundary)
%   permeability: double(1, nboundary)
%
%   flux_matrix: double(nnode, nnode)


% Number of boundaries
nboundary = length(facets);

% Initialize output matrix
flux_matrix = sparse(length(points), length(points));

% Add block in matrix for each touching boundary
for iboundary = 1:nboundary
    
    % Check that the boundary is permeable
    if permeability(iboundary) > 1e-16
        
        % Extract boundary facets
        boundary = facets{iboundary}';
        
        % Only proceed if the boundary touches the current compartment
        if ~isempty(boundary)
            
            % Identify nodes in boundary
            inds = unique(boundary);

            % Set weigths to boundary permeability for boundary nodes, else 0
            coeffs = zeros(max(inds), 1);
            coeffs(inds) = permeability(iboundary);

            % Add block to flux matrix (the blocks do not overlap)
            flux_matrix = flux_matrix + flux_matrixP1_3D(boundary, points', coeffs);
        end
    end
end
