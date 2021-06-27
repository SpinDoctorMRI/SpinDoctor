function plot_geometry_info(setup, femesh)
%PLOT_GEOMETRY_INFO Plot information of the geometry.
%
%   setup: struct
%   volumes: [1 x ncompartment]
%   surface_areas
%
%   1 figure with title of "Connections boundary-compartment"
%   1 figure with title of "volumes"
%   1 figure with title of "Surface Area"

compartments = setup.pde.compartments;

cmpts_in = compartments == "in";
cmpts_out = compartments == "out";
cmpts_ecs = compartments == "ecs";

% Get volume and surface area quantities from mesh
[volumes, surface_areas] = get_vol_sa(femesh);
boundary_markers = ~cellfun(@isempty, femesh.facets);

ncompartment = length(volumes);
nboundary = size(boundary_markers, 2);

figure;
hold all;

cmpts_boundaries_mat_plot = zeros(size(boundary_markers));
cmpts_boundaries_mat_plot(cmpts_in, :) = boundary_markers(cmpts_in, :);
spy(cmpts_boundaries_mat_plot, "kd", 10);

cmpts_boundaries_mat_plot(:, :) = 0;
cmpts_boundaries_mat_plot(cmpts_out, :) = boundary_markers(cmpts_out, :);
spy(cmpts_boundaries_mat_plot, "ko", 10);

cmpts_boundaries_mat_plot(:, :) = 0;
cmpts_boundaries_mat_plot(cmpts_ecs, :) = boundary_markers(cmpts_ecs, :);
spy(cmpts_boundaries_mat_plot, "ks", 10);


xlabel("Boundary");
ylabel("Compartment");
title("Connections boundary-compartment");
set(gca, "xtick", 1:nboundary);
set(gca, "ytick", 1:ncompartment);
grid on

xticklabel = num2cell(1:ncompartment);
xticklabel{end + 1} = "Total";

figure;
hold on;
bar(1:ncompartment, volumes, "b");
bar(ncompartment + 1, sum(volumes), "r");
title("Volumes");
set(gca, "xtick", 1:ncompartment+1);
set(gca, "xticklabel", xticklabel);
xlabel("Compartment");
grid on

figure;
hold on;
bar(1:ncompartment, surface_areas, "b");
bar(ncompartment + 1, sum(surface_areas), "r");
title("Surface Area");
set(gca, "xtick", 1:ncompartment+1);
set(gca, "xticklabel", xticklabel);
xlabel("Compartment");
grid on
