function femesh = split_mesh(femesh_all)
%SPLIT_MESH Split mesh into compartments.
%
%   femesh_all: struct with fields
%       points_all: [3 x npoint]
%       facets_all: [3 x nfacet]
%       elements_all: [4 x nelement]
%       facetmarkers: [1 x nfacet]
%       elementmarkers: [1 x nelement]
%
%   femesh: struct with fields
%       ncompartment: [1 x 1]
%         nboundary: [1 x 1]
%       points: {1 x ncompartment}[3 x npoint]
%       facets: {ncompartment x nboundary}[3 x nfacet]
%       elements: {1 x ncompartment}[4 x nelement]
%       point_map: {1 x ncompartment}[npoint x 1]


% Extract global mesh
points_all = femesh_all.points;
facets_all = femesh_all.facets;
elements_all = femesh_all.elements;
facetmarkers = femesh_all.facetmarkers;
elementmarkers = femesh_all.elementmarkers;

% Identify compartments and boundaries
compartments = unique(elementmarkers);
boundaries = unique(facetmarkers);
ncompartment = length(compartments);
nboundary = length(boundaries);

% Split points and elements into compartment (with ghost points)
fprintf("Separating FE mesh into %d compartments\n", ncompartment);
elements = cell(1, ncompartment);
point_map = cell(1, ncompartment);
points = cell(1, ncompartment);
for icmpt = 1:ncompartment
    elements{icmpt} = elements_all(:, elementmarkers == compartments(icmpt));
    point_map{icmpt} = unique(elements{icmpt});
    points{icmpt} = points_all(:, point_map{icmpt});
end

% Split facets into boundaries
fprintf("Separating FE mesh with %d boundaries\n", nboundary);
boundary_facets = cell(1, nboundary);
for iboundary = 1:nboundary
    boundary_facets{iboundary} = facets_all(:, facetmarkers == boundaries(iboundary));
end

% Make sure elements facets and boundary points refer to the points in the new
% numbering system
facets = cell(ncompartment, nboundary);

for icmpt = 1:ncompartment
    % Renumber nodes
    oldcode = point_map{icmpt};
    newcode = 1:length(point_map{icmpt});
    assert(numel(newcode) == numel(oldcode),  ...
        "newcode and oldcode must have the same number of elements");
    [toreplace, bywhat] = ismember(elements{icmpt}, oldcode);
    elements{icmpt}(toreplace) = newcode(bywhat(toreplace));

    for iboundary = 1:nboundary
        % Check that the boundary lies on the compartment in its entirety
        boundary_on_compartment = all(ismember(boundary_facets{iboundary}, point_map{icmpt}), "all");
        if boundary_on_compartment
            % Associate boundary to compartment in new local numbering system

            % Create list of facets in boundary in new local numbering system
            facets{icmpt, iboundary} = boundary_facets{iboundary};
            oldcode = point_map{icmpt};
            newcode = 1:length(point_map{icmpt});
            assert(numel(newcode) == numel(oldcode), ...
                "newcode and oldecode must have the same number of elements");
            [toreplace, bywhat] = ismember(facets{icmpt, iboundary}, oldcode);
            facets{icmpt, iboundary}(toreplace) = newcode(bywhat(toreplace));
        end
    end
end

% Checking the mesh quality
hmax = 0;
% aspect_ratio_lim = 0.05; % the worst -> [0, 1] <- the best
for icmpt = 1:ncompartment
    qmesh = mesh_quality(points{icmpt}, elements{icmpt});
    hmax = max(max(qmesh.hout), hmax);
    fprintf("  Compartment %d of %d – FE mesh with minimum aspect ratio of %.1e\n", ...
        icmpt, ncompartment, qmesh.quality(1));
    % if qmesh{icmpt}.quality1(1) < aspect_ratio_lim
    %     fprintf("  Compartment %d of %d – FE mesh with minimum aspect ratio of %.1e\n", ...
    %         icmpt, ncompartment, qmesh.quality1(1));
    % end
end

% Create output structure
femesh.ncompartment = ncompartment;
femesh.nboundary = nboundary;
femesh.points = points;
femesh.facets = facets;
femesh.elements = elements;
femesh.point_map = point_map;
