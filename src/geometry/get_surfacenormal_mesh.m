function [areas, facet_centers, normals] = get_surfacenormal_mesh(points, elements, facets)
%GET_SURFACENORMAL_MESH Start Computing normals FacN from given facets facets.


[~, areas, facet_centers, normals] = get_surface_mesh(points, facets);

% Find elements to which the facets belong
neumann2element = match_boundary_facets_with_elements(facets', elements');
elements_match = elements(:, neumann2element);
nelement_match = size(elements_match, 2);

% Identify centers of matching elements
x = reshape(points(:, elements_match), 3, 4, nelement_match);
element_centers = squeeze(mean(x, 2));

% Create outward directed vectors (from cell center to outer facet)
outwards_vectors = facet_centers - element_centers;

% Determine orientation of normals (+ for out, - for in)
orientations = sign(dot(outwards_vectors, normals));

% Orient all normals outwards
normals = orientations .* normals;

% figure;
% hold on;
% h = quiver3(facet_centers(1, :), facet_centers(2, :), facet_centers(3, :), ...
%     normals(1, :), normals(2, :), normals(3, :), 0, "k");
% h = plot3(element_centers(1, :), element_centers(2, :), element_centers(3, :), "ro");
% view(3);
% axis equal;
