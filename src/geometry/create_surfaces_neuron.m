function surfaces = create_surfaces_neuron(filename, setup)
%CREATE_SURFACES_NEURON Create neuron surface mesh.
%   A neuron surface mesh is loaded or create and loaded. An ECS can be added.
%
%   filename: string
%   setup: struct
%
%   surfaces: struct with fields
%       points: [3 x npoint]
%       facets: [3 x nfacet]
%       facetmarkers: [1 x nfacet]
%       regions: [3 x ncompartment]


% Number of points to discretize space for creating tight wrap ECS
ndiscretize = 100;

ecs_shape = setup.geometry.ecs_shape;
include_ecs = ecs_shape ~= "no_ecs";
if include_ecs
    ecs_ratio = setup.geometry.ecs_ratio;
else
    % the default ecs_ratio is used to find a point inside neuron
    ecs_ratio = 0.1;
end

if endsWith(filename, ".1")
    femesh = read_tetgen(filename);
    points = femesh.points;
    elements = femesh.elements;
    facets = femesh.facets;
else
    if ~isfile(filename + "_elements.txt") || ~isfile(filename + "_nodes.txt")
        gmsh_to_fem_try(filename);
    end
    disp("Reading from neuron FE mesh from " + filename);

    elements = load(filename + "_elements.txt");
    points = load(filename + "_nodes.txt");
    points = points';
    elements = elements';
    
    [elems2facets, facets2points] = get_faces(elements');
    facets2elems = entryInWhichRows(elems2facets);
    bindex = facets2elems(:, 2) == 0;  % index of boundary facets - belong to one element only
    bb = facets2points(bindex, :);     % points of boundary facets
    aa = sort(bb, 2);
    facets = unique(aa, "rows");
    facets = facets';
end

facetmarkers = ones(1, size(facets, 2));
[regions, regions_ecs] = find_regions(points, elements, facets, ecs_ratio);

pmin = min(points, [], 2);
pmax = max(points, [], 2);

if ecs_shape ~= "no_ecs"
    ecs_gap = ecs_ratio * 10;

    switch ecs_shape
        case "box"
            pmin = pmin - ecs_gap;
            pmax = pmax + ecs_gap;
            points_ecs = [
                pmin(1) pmax(1) pmax(1) pmin(1) pmin(1) pmax(1) pmax(1) pmin(1)
                pmin(2) pmin(2) pmax(2) pmax(2) pmin(2) pmin(2) pmax(2) pmax(2)
                pmin(3) pmin(3) pmin(3) pmin(3) pmax(3) pmax(3) pmax(3) pmax(3)
            ];
            facets_ecs = [
                1 1 1 1 2 2 3 3 4 4 5 5
                2 3 2 6 3 7 4 8 1 5 6 7
                3 4 6 5 7 6 8 7 5 8 7 8
            ];
        case "convex_hull"
            DT = delaunayTriangulation(points');
            [T, Xb] = freeBoundary(DT);
            TR = triangulation(T, Xb);
            x1 = TR.Points;
            x2 = incenter(TR);
            n1 = TR.vertexNormal;
            n2 = TR.faceNormal;
            y1 = x1 + ecs_gap * n1;
            y2 = x2 + ecs_gap * n2;
            y = [y1; y2];
            DT = delaunayTriangulation(y);

            pmin = pmin - 10 * ecs_gap;
            pmax = pmax + 10 * ecs_gap;
            p = pmin + (pmax - pmin) .* linspace(0, 1, ndiscretize);
            [X, Y, Z] = meshgrid(p(1, :), p(2, :), p(3, :));
            p = [X(:) Y(:) Z(:)];
            inds = DT.pointLocation(p);
            markers = ones(size(inds));
            markers(isnan(inds)) = -0.2;
            V = reshape(markers, ndiscretize, ndiscretize, ndiscretize);
            FV = isosurface(X, Y, Z, V, 0);
            points_ecs = FV.vertices';

            shp_ecs = alphaShape(points_ecs');
            shp_ecs.HoleThreshold = prod(pmax - pmin);
            shp_ecs.Alpha = 2.5 * shp_ecs.criticalAlpha("one-region");
            [bf , P] = boundaryFacets(shp_ecs);
            points_ecs = P';
            facets_ecs = bf';
        case "tight_wrap"
            shp = alphaShape(points', 'HoleThreshold', (4/3)*pi*20^3);
            shp.Alpha = 1.2*shp.criticalAlpha("one-region");
            
            % normal vector method to reconstruct ECS surface mesh
            [tri, xyz] = boundaryFacets(shp);
            TR = triangulation(tri,xyz);
            V = vertexNormal(TR);
            P = incenter(TR);
            F = faceNormal(TR);
            
            points_ecs = [(xyz + 0.5*ecs_gap * V)',...
                          (xyz + ecs_gap * V)',...
                          (P + 0.5*ecs_gap * F)',...
                          (P + ecs_gap * F)'];
            shp2 = alphaShape([points,points_ecs]', 'HoleThreshold', (4/3)*pi*20^3);
            shp2.Alpha = 1.5*shp2.criticalAlpha("one-region");
            [tri, xyz] = boundaryFacets(shp2);
            points_ecs = xyz';
            facets_ecs = tri';
            
%             % isosurface method to reconstruct ECS surface mesh
%             pmin = pmin - 2 * ecs_gap;
%             pmax = pmax + 2 * ecs_gap;
%             p = pmin + (pmax - pmin) .* linspace(0, 1, ndiscretize);
%             [X, Y, Z] = meshgrid(p(1, :), p(2, :), p(3, :));
%             p = [X(:) Y(:) Z(:)];
%             inside = shp.inShape(p);
%             [~, dist] = shp.nearestNeighbor(p);
%             dist(inside) = -dist(inside);
%             V = reshape(dist, ndiscretize, ndiscretize, ndiscretize);
%             FV = isosurface(X, Y, Z, V, ecs_gap);
%             points_ecs = FV.vertices';
%             facets_ecs = FV.faces';
    end
    
    facetmarkers_ecs = 2 * ones(1, size(facets_ecs, 2));
    
    npoint_out = size(points, 2);
    points = [points points_ecs];
    facets = [facets facets_ecs+npoint_out];
    facetmarkers = [facetmarkers facetmarkers_ecs];
    regions = [regions regions_ecs];
end


% Create output structure
surfaces.points = points;
surfaces.facets = facets;
surfaces.facetmarkers = facetmarkers;
surfaces.regions = regions;
