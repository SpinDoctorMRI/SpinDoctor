function [total_area, areas, centers, normals] = get_surface_mesh(points, facets)
%GET_SURFACE_MESH Compute surface areas, centers and normals for each facet.
%
%   points: double(3, npoint)
%   facets: int(3, nfacet)
%
%   total_area: double
%   areas: double(1, nfacet)
%   centers: double(3, nfacet)
%   normals: double(3, nfacet)


% Number of facets
nfacet = size(facets, 2);

% Facets
x = reshape(points(:, facets), 3, 3, nfacet);

% Facet centers
centers = squeeze(mean(x, 2));

% Facet normals
normals = squeeze(cross(x(:, 1, :) - x(:, 2, :), x(:, 3, :) - x(:, 2, :)));

% Facet areas
areas = 1 / 2 * vecnorm(normals);
total_area = sum(areas);

% Normalize normals
normals = normals ./ vecnorm(normals);
