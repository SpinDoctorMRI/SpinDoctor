function surfaces = create_surfaces_sphere(cells, params_cells)
%CREATE_SURFACES_SPHERE Create sphere surface mesh.
%   Surfaces may include in, out, and ecs compartments.
%
%   cells: struct with fields
%       centers: [3 x ncell]
%       radii: [1 x ncell]
%   params_cells: struct
%
%   surfaces: struct with fields
%       points: [3 x npoint]
%       facets: [3 x nfacet]
%       facetmarkers: [1 x nfacet]
%       regions: [3 x ncompartment]

% Number of points for an average sphere
nmin = 50;
% nref = 200;
nref = 300;
% nref_ecs = 100;

% Number of points to discretize space for creating tight wrap ECS
ndiscretize = 50;

% Extract parameters
centers = cells.centers;
radii = cells.radii;
ncell = params_cells.ncell;
rmin = params_cells.rmin;
rmax = params_cells.rmax;
include_in = params_cells.include_in;
in_ratio = params_cells.in_ratio;
ecs_shape = params_cells.ecs_shape;
ecs_ratio = params_cells.ecs_ratio;

include_ecs = ecs_shape ~= "no_ecs";
rmean = (rmin + rmax) / 2;

create_npoint = @(r) max(nmin, round(nref * (r / rmean).^2));

if include_in
    npoint_in = create_npoint(in_ratio .* radii);
    points_in = cell(1, ncell);
    facets_in = cell(1, ncell);
    for icell = 1:ncell
        p = centers(:, icell) + in_ratio * radii(icell) .* fibonacci(npoint_in(icell));
        points_in{icell} = p;
        facets_in{icell} = convhull(p')';
    end
    regions_in = centers;
    nfacet_in = cellfun(@(x) size(x, 2), facets_in);
else
    npoint_in = zeros(1, 0);
    nfacet_in = zeros(1, 0);
    points_in = {};
    facets_in = {};
    regions_in = zeros(3, 0);
end

npoint_out = create_npoint(radii);
points_out = cell(1, ncell);
facets_out = cell(1, ncell);
for icell = 1:ncell
    p = centers(:, icell) + radii(icell) .* fibonacci(npoint_out(icell));
    points_out{icell} = p;
    facets_out{icell} = convhull(p')';
end
regions_out = centers + [1; 0; 0] * (include_in * in_ratio + 1) / 2 * radii;
nfacet_out = cellfun(@(x) size(x, 2), facets_out);

if include_ecs
    radii_ecs = radii + ecs_ratio * rmean;
    npoint_ecs = create_npoint(radii_ecs);
    % npoint_ecs = repelem(nref_ecs, 1, ncell); 
    points_ecs = cell(1, ncell);
    for icell = 1:ncell
        p = centers(:, icell) + radii_ecs(icell) .* fibonacci(npoint_ecs(icell));
        points_ecs{icell} = p;
    end
    points_ecs = [points_ecs{:}];
    [~, ind] = sort(centers(1, :));
    regions_ecs = centers(:, ind(1)) - [1; 0; 0] * (radii(ind(1)) + ecs_ratio / 2 * rmean);
    switch ecs_shape
        case "box"
            pmin = min(points_ecs, [], 2);
            pmax = max(points_ecs, [], 2);
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
%             facets_ecs = convhull(points_ecs')';
%             inds = unique(facets_ecs);
%             points_ecs = points_ecs(:, inds);
%             tmp = facets_ecs == shiftdim(inds, -2);
%             [i, j, k] = ind2sub(size(tmp), find(tmp));
%             facets_ecs(sub2ind(size(facets_ecs), i, j)) = k;

            DT = delaunayTriangulation(points_ecs');
            
            % [T, Xb] = freeBoundary(DT);
            % TR = triangulation(T, Xb);
            % x1 = TR.Points;
            % x2 = TR.incenter;
            % n1 = TR.vertexNormal;
            % n2 = TR.faceNormal;
            % y1 = x1 + ecs_ratio * rmean * n1;
            % y2 = x2 + ecs_ratio * rmean * n2;
            % y = [y1; y2];
            % k = convhull(y);
            % points_ecs = y';
            % facets_ecs = k';
            
            pmin = min(points_ecs, [], 2) - 3 * ecs_ratio * rmax;
            pmax = max(points_ecs, [], 2) + 3 * ecs_ratio * rmax;
            p = pmin + (pmax - pmin) .* linspace(0, 1, ndiscretize);
            [X, Y, Z] = meshgrid(p(1, :), p(2, :), p(3, :));
            p = [X(:) Y(:) Z(:)];
            inds = DT.pointLocation(p);
            markers = ones(size(inds));
            markers(isnan(inds)) = -0.2;
            V = reshape(markers, ndiscretize, ndiscretize, ndiscretize);
            V = smooth3(V, "gaussian", [9, 9, 9]);
            FV = isosurface(X, Y, Z, V, 0);
            points_ecs = FV.vertices';
            facets_ecs = FV.faces';
            
            shp_ecs = alphaShape(points_ecs');
            shp_ecs.HoleThreshold = prod(pmax - pmin);
            shp_ecs.Alpha = 2 * shp_ecs.criticalAlpha("one-region");
            % shp_ecs.Alpha = max(pmax - pmin);
            [bf , P] = boundaryFacets(shp_ecs);
            points_ecs = P';
            facets_ecs = bf';
            
        case "tight_wrap"
            % shp_ecs = alphaShape(points_ecs');
            % shp_ecs.HoleThreshold = 2 * 4 * pi / 3 * max(radii_ecs).^2;
            % shp_ecs.Alpha = 1.05 * shp_ecs.criticalAlpha("one-region");
            % [bf , P] = boundaryFacets(shp_ecs);
            % points_ecs = P';
            % facets_ecs = bf';
            
            pmin = min(points_ecs, [], 2) - ecs_ratio * rmean;
            pmax = max(points_ecs, [], 2) + ecs_ratio * rmean;
            x = linspace(pmin(1, :), pmax(1, :), ndiscretize);
            y = linspace(pmin(1, :), pmax(1, :), ndiscretize);
            z = linspace(pmin(1, :), pmax(1, :), ndiscretize);
            [X, Y, Z] = meshgrid(x, y, z);
            p = [X(:) Y(:) Z(:)]';
            distmats = vecnorm(p - permute(centers, [1, 3, 2]));
            distmats = distmats - shiftdim(radii, -1);
            distmats = max(0, distmats);
            distmats = min(distmats, [], 3);
            V = reshape(distmats, ndiscretize, ndiscretize, ndiscretize);
            FV = isosurface(X, Y, Z, V, ecs_ratio * rmax);
            points_ecs = FV.vertices';
            facets_ecs = FV.faces';
    end
    npoint_ecs = size(points_ecs, 2);
    nfacet_ecs = size(facets_ecs, 2);
else
    npoint_ecs = zeros(1, 0);
    nfacet_ecs = zeros(1, 0);
    points_ecs = zeros(3, 0);
    facets_ecs = zeros(3, 0);
    regions_ecs = zeros(3, 0);
end

% Merge triangulations
npoint = cumsum([0, npoint_in, npoint_out, npoint_ecs]);
nfacet = cumsum([0, nfacet_in, nfacet_out, nfacet_ecs]);
points = [points_in{:}, points_out{:}, points_ecs];
facets = [facets_in{:}, facets_out{:}, facets_ecs];
facetmarkers = zeros(1, size(facets, 2));
for i = 1:length(npoint)-1
    inds = nfacet(i)+1:nfacet(i+1);
    facets(:, inds) = facets(:, inds) + npoint(i);
    facetmarkers(inds) = i;
end
regions = [regions_in regions_out regions_ecs];


% Create output structure
surfaces.points = points;
surfaces.facets = facets;
surfaces.facetmarkers = facetmarkers;
surfaces.regions = regions;
