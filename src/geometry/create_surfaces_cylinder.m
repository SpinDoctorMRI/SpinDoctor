function surfaces = create_surfaces_cylinder(cells, setup)
%CREATE_SURFACES_CYLINDER Create cylinder surface mesh.
%   Surfaces may include in, out, and ecs compartments.
%
%   cells: struct with fields
%       centers: [2 x ncell]
%       radii: [1 x ncell]
%   setup: struct
%
%   surfaces: struct with fields
%       points: [3 x npoint]
%       facets: [3 x nfacet]
%       facetmarkers: [1 x nfacet]
%       regions: [3 x ncompartment]


% Design parameters
nside = 30;
nside_min = 12;

% Extract parameters
centers = cells.centers;
radii = cells.radii;
ncell = setup.geometry.ncell;
rmin = setup.geometry.rmin;
rmax = setup.geometry.rmax;
height = setup.geometry.height;
include_in = setup.geometry.include_in;
in_ratio = setup.geometry.in_ratio;
ecs_shape = setup.geometry.ecs_shape;
ecs_ratio = setup.geometry.ecs_ratio;

include_ecs = ecs_shape ~= "no_ecs";
nboundary = (2 * include_in + 1 + include_ecs) * ncell + include_ecs;

rmean = (rmin + rmax) / 2;

create_nside = @(r) max(nside_min, round(nside * r / rmean));
% create_nside = @(r) nside;
create_angles = @(n) 2 * pi * (0:n-1) / n;
create_circle = @(center, r, angles) center + r .* [cos(angles); sin(angles)];
create_edges = @(n) [1:n-1 n; 2:n 1];


if include_in
    radii_in = in_ratio * radii;
    nside_in = create_nside(radii_in);
    circles_in = cell(1, ncell);
    circle_edges_in = cell(1, ncell);
    for icell = 1:ncell
        circles_in{icell} = create_circle(centers(:, icell), radii_in(icell), ...
            create_angles(nside_in(icell)));
        circle_edges_in{icell} = create_edges(nside_in(icell));
    end
    points_in = [circles_in{:}];
    edges_in = [circle_edges_in{:}];
    nold = 0;
    for icell = 1:ncell
        n = nold + nside_in(icell);
        edges_in(:, nold+1:n) = edges_in(:, nold+1:n) + nold;
        nold = n;
    end
    regions_in = centers;
else
    nside_in = zeros(1, 0);
    points_in = zeros(2, 0);
    edges_in = zeros(2, 0);
    regions_in = zeros(2, 0);
end

nside_out = create_nside(radii);
circles_out = cell(1, ncell);
circle_edges_out = cell(1, ncell);
for icell = 1:ncell
    circles_out{icell} = create_circle(centers(:, icell), radii(icell), ...
        create_angles(nside_out(icell)));
    circle_edges_out{icell} = create_edges(nside_out(icell));
end
points_out = [circles_out{:}];
edges_out = [circle_edges_out{:}];
nold = 0;
for icell = 1:ncell
    n = nold + nside_out(icell);
    edges_out(:, nold+1:n) = edges_out(:, nold+1:n) + nold;
    nold = n;
end
regions_out = centers + [1; 0] * (include_in * in_ratio + 1) / 2 * radii;

if include_ecs
    radii_ecs = radii + ecs_ratio * rmean;
    nside_ecs = create_nside(radii_ecs);
    circles_ecs = cell(1, ncell);
    for icell = 1:ncell
        circles_ecs{icell} = create_circle(centers(:, icell), radii_ecs(icell), ...
            create_angles(nside_ecs(icell)));
    end
    points_ecs = [circles_ecs{:}];
    [~, ind] = sort(centers(1, :));
    regions_ecs = centers(:, ind(1)) - [1; 0] * (radii(ind(1)) + ecs_ratio / 2 * rmean);

    switch ecs_shape
        case "box"
            pmin = min(points_ecs, [], 2);
            pmax = max(points_ecs, [], 2);
            points_ecs = [
                pmin(1) pmax(1) pmax(1) pmin(1)
                pmin(2) pmin(2) pmax(2) pmax(2)
            ];
            edges_ecs = [
                1 2 3 4
                2 3 4 1
            ];
        case "convex_hull"
            inds = convhull(points_ecs')';
            points_ecs = points_ecs(:, inds(1:end-1));
            edges_ecs = [1:length(inds)-1; 2:length(inds)-1 1];
        case "tight_wrap"
            shp_ecs = alphaShape(points_ecs');
            shp_ecs.HoleThreshold = 2 * pi * ((1 + ecs_ratio) * max(radii)).^2;
            shp_ecs.Alpha = 2.5 * shp_ecs.criticalAlpha("one-region");
            [bf, P] = boundaryFacets(shp_ecs);
            points_ecs = P';
            edges_ecs = bf';
    end
    nedge_ecs = size(edges_ecs, 2);
else
    points_ecs = zeros(2, 0);
    edges_ecs = zeros(2, 0);
    regions_ecs = zeros(2, 0);
    nedge_ecs = [];
end


npoint_in = size(points_in, 2);
npoint_out = size(points_out, 2);
npoint_ecs = size(points_ecs, 2);



points = [points_in points_out points_ecs];
edges = [edges_in edges_out+npoint_in edges_ecs+npoint_in+npoint_out];

npoint = size(points, 2);

if ecs_shape == "tight_wrap"
    shp = alphaShape(points');
    shp.Alpha = shp_ecs.Alpha;
    shp.HoleThreshold = 2 * pi * (2 * (1 + ecs_ratio) * max(radii)).^2;
    facets = shp.alphaTriangulation';
    [~, P] = boundaryFacets(shp);
    if size(P, 1) ~= size(points_ecs, 2)
        plot(shp_ecs);
        hold on;
        plot(shp, "FaceColor", "r", "FaceAlpha", 0.5, "EdgeColor", "b")
        title("Error in green area.")
        error("Problem in cylinder ground mesh tight wrap generation. " ...
            + "Try deleting geometry files and rerun geometry generation.");
    end
else
    DT = delaunayTriangulation(points', edges');
    facets = DT.ConnectivityList';
end

nfacet = size(facets, 2);


boundary_bounds = cumsum([0 nside_in nside_out nedge_ecs]);

ncell_in = include_in * ncell;

find_boundary = @(node) find(boundary_bounds(1:end-1) + 1 <= node & node <= boundary_bounds(2:end), 1);

% Create boundary numbers
facetmarkers = zeros(1, nfacet);
for ifacet = 1:nfacet
    b = arrayfun(find_boundary, facets(:, ifacet));
    b_in = b <= ncell_in;
    b_out = ncell_in < b & b <= ncell_in + ncell;
    if all(b_in)
        % Inside IN compartment
        facetmarkers(ifacet) = ncell_in + include_ecs * ncell + b(1);
    elseif any(b_in)
        % Triangle touches IN compartment, but is not inside. It is thus the corresponding OUT compartment.
        tmp = b(~b_in);
        facetmarkers(ifacet) = ncell_in + include_ecs * ncell + tmp(1);
    elseif all(b_out) && b(1) == b(2) && b(1) == b(3)
        % Triangle lies fully within the same OUT comparment
        facetmarkers(ifacet) = ncell_in + include_ecs * ncell + b(1);
    else
        % Triangle lies in the ECS
        if include_ecs
            facetmarkers(ifacet) = nboundary;
        else
            facetmarkers(ifacet) = NaN;
        end
    end
end

if ~include_ecs
    % remove ECS facets, if no_ecs
    facets = facets(:, ~isnan(facetmarkers));
    facetmarkers = facetmarkers(~isnan(facetmarkers));
    nfacet = size(facets, 2);
end

% Copy points, edges and markers two the two 3D planes top and bottom
z = repelem(height / 2, 1, npoint);
points = [points points; -z z];

create_sidefacets = @(le, ri) [
   (le:ri-1)        ri        (le+1:ri)+npoint le+npoint
   (le:ri-1)+npoint ri+npoint (le:ri-1)+npoint ri+npoint
   (le+1:ri)        le        (le+1:ri)        le
];

sidefacets = arrayfun(create_sidefacets, boundary_bounds(1:end-1) + 1, ...
    boundary_bounds(2:end), "UniformOutput", false);
sidefacets = [sidefacets{:}];

sidefacetmarkers_in = 1:ncell_in;
sidefacetmarkers_in = arrayfun(@(m, n) repmat([m m], 1, n), sidefacetmarkers_in, ...
    nside_in, "UniformOutput", false);
sidefacetmarkers_in = [sidefacetmarkers_in{:}];

sidefacetmarkers_out = (1 + ~include_ecs)*ncell_in+1:(1 + ~include_ecs)*ncell_in+ncell;
sidefacetmarkers_out = arrayfun(@(m, n) repmat([m m], 1, n), sidefacetmarkers_out, ...
    nside_out, "UniformOutput", false);
sidefacetmarkers_out = [sidefacetmarkers_out{:}];

sidefacetmarkers_ecs = repelem(nboundary, 1, 2*npoint_ecs);

facets = [facets facets+npoint sidefacets];
facetmarkers = [facetmarkers facetmarkers sidefacetmarkers_in sidefacetmarkers_out sidefacetmarkers_ecs];

regions = [regions_in regions_out regions_ecs];
regions = [regions; zeros(1, size(regions, 2))];


% Create output structure
surfaces.points = points;
surfaces.facets = facets;
surfaces.facetmarkers = facetmarkers;
surfaces.regions = regions;
