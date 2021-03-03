function plot_surface_triangulation(surfaces)
%PLOT_SURFACE_TRIANGULATION Plot surface triangulation.
%


p = surfaces.points;
facets = surfaces.facets;
% facetmarkers = surfaces.facetmarkers;
regions = surfaces.regions;

figure;
hold on
trisurf(facets', p(1, :), p(2, :), p(3, :), "facecolor", "interp", "facealpha", 0.6)
plot3(regions(1, :), regions(2, :), regions(3, :), "*");
axis equal;
view(3);
title("Surface triangulation of canonical configuration");
