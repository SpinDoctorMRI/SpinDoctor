function plot_femesh(femesh, compartments)
%PLOT_FEMESH Plot the finite element meshes of inner, outer, and ECS compartments.
%
%   femesh: struct
%   compartments


ncompartment = length(compartments);
include_in = any(compartments == "in");
include_ecs = any(compartments == "ecs");

cmpts = "out";
figs.out = figure;
hold on
if include_in
    cmpts = ["in" cmpts];
    figs.in = figure;
    hold on
end
if include_ecs
    cmpts = [cmpts "ecs"];
    figs.ecs = figure;
    hold on
end

% Determine limits of domain
pmin = min([femesh.points{:}], [], 2);
pmax = max([femesh.points{:}], [], 2);
axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];

for icmpt = 1:ncompartment
    facets = [femesh.facets{icmpt, :}];
    points = femesh.points{icmpt};
    
    figure(figs.(compartments(icmpt)));
    h = trisurf(facets', points(1, :), points(2, :), points(3, :));
    set(h, "facealpha", 0.7);
    set(h, "facecolor", "interp");
    % set(h, "LineStyle", "none");
    set(h, "LineWidth", 0.03);
    center = mean(points, 2);
    center(3) = max(points(3, :)) + 1;
    text(center(1), center(2), center(3), num2str(icmpt), "color", "r", "fontsize", 20);
end

for cmpt = cmpts
    figure(figs.(cmpt));
    xlabel("x");
    ylabel("y");
    zlabel("z");
    grid on;
    view(3);
    axis equal;
    axis(axis_vec);
    inds = find(compartments == cmpt);
    cmpt_str = sprintf(join(repmat("%d", 1, length(inds))), inds);
    title(sprintf("Finite element mesh, %s compartment: [%s]", cmpt, cmpt_str));
end
