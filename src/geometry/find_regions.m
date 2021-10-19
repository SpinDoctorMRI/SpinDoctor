function [regions, regions_ecs] = find_regions(p, e, t, ecs_ratio)
    [~, facet_centers, normals] = get_surfacenormal_mesh(p, e, t);
    ind = randi(size(facet_centers, 2));
    regions = facet_centers(:, ind) - ecs_ratio*normals(:, ind)/2;
    regions_ecs = facet_centers(:, ind) + ecs_ratio*normals(:, ind)/2;
end