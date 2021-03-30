function flux_matrices = assemble_flux_matrix(points, facets)
%ASSEMBLE_FLUX_MATRIX Assemble flux matrix in compartment
%
%   points: {1 x ncompartment}[3 x npoint]
%   facets: {ncompartment x nboundary}[[3 x nfacet]
%
%   flux_matrix: {ncompartment x nboundary}[npoint x npoint]


% Number of boundaries
[ncompartment, nboundary] = size(facets);

% Initialize output matrices
flux_matrices = cell(ncompartment, nboundary);

for icmpt = 1:ncompartment
    for iboundary = 1:nboundary
        % Extract boundary facets
        boundary = facets{icmpt, iboundary};
        
        % Only proceed if the boundary touches the current compartment
        if ~isempty(boundary)
            % Assemble flux matrix compartment-boundary
            flux_matrices{icmpt, iboundary} = flux_matrixP1_3D(boundary', points{icmpt}');
        end
    end
end
