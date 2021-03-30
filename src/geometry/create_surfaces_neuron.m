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
ecs_ratio = setup.geometry.ecs_ratio;

if ~isfile(filename + "_elements.txt") || ~isfile(filename + "_nodes.txt")
    gmsh_to_fem_try(filename);
end
disp("Reading from neuron FE mesh from " + filename);

elements = load(filename + "_elements.txt");
points = load(filename + "_nodes.txt");
points = points';
elements = elements';

[total_volume, volumes, centers] = get_volume_mesh(points, elements);
centermass = centers * volumes' / total_volume;

[elems2facets, facets2points] = get_faces(elements');
facets2elems = entryInWhichRows(elems2facets);
bindex = facets2elems(:, 2) == 0; % index of boundary facets - belong to one element only
bb = facets2points(bindex, :);     % points of boundary facets
aa = sort(bb, 2);
facets = unique(aa, "rows");
facets = facets';
facetmarkers = ones(1, size(facets, 2));
regions = centermass;

pmin = min(points, [], 2);
pmax = max(points, [], 2);

if ecs_shape ~= "no_ecs"
    
%     ecs_gap = ecs_ratio * min(pmax - pmin);
    ecs_gap = ecs_ratio * 10;

    switch ecs_shape
        case "box"
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
            V = smooth3(V, "gaussian");%, [9, 9, 9]);
            FV = isosurface(X, Y, Z, V, 0);
            points_ecs = FV.vertices';
            % facets_ecs = FV.faces';
            
            shp_ecs = alphaShape(points_ecs');
            shp_ecs.HoleThreshold = prod(pmax - pmin);
            shp_ecs.Alpha = 2 * shp_ecs.criticalAlpha("one-region");
            % shp_ecs.Alpha = max(pmax - pmin);
            [bf , P] = boundaryFacets(shp_ecs);
            points_ecs = P';
            facets_ecs = bf';
        case "tight_wrap"
            shp = alphaShape(points');
            shp.Alpha = 1.2 * shp.criticalAlpha("one-region");
            
            pmin = pmin - 2 * ecs_gap;
            pmax = pmax + 2 * ecs_gap;
            p = pmin + (pmax - pmin) .* linspace(0, 1, ndiscretize);
            [X, Y, Z] = meshgrid(p(1, :), p(2, :), p(3, :));
            p = [X(:) Y(:) Z(:)];
            
            inside = shp.inShape(p);
            [~, dist] = shp.nearestNeighbor(p);
            dist(inside) = -dist(inside);

            V = reshape(dist, ndiscretize, ndiscretize, ndiscretize);
%             V = smooth3(V, "gaussian");%, [9, 9, 9]);
            FV = isosurface(X, Y, Z, V, ecs_gap);
            points_ecs = FV.vertices';
            facets_ecs = FV.faces';
    end
    
    facetmarkers_ecs = 2 * ones(1, size(facets_ecs, 2));
    [~, ind] = sort(points(1, :));
    ind = ind(1);
    regions_ecs = points(:, ind) - ecs_gap / 10 * [1; 0; 0];
    
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
