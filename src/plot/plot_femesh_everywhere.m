function plot_femesh_everywhere(femesh, subtitle_str)
%PLOT_FEMESH_EVERYWHERE Plot the finite element mesh of all compartments.
%
%   femesh: struct


% Number of compartments
ncompartment = femesh.ncompartment;

% Determine limits of domain
pmin = min([femesh.points{:}], [], 2);
pmax = max([femesh.points{:}], [], 2);
axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];

figure;
hold on;
for icmpt = [ncompartment 1:ncompartment-1] % Plot ECS first
    facets = [femesh.facets{icmpt, :}];
    points = femesh.points{icmpt};
    h = trisurf(facets', points(1, :), points(2, :), points(3, :));
    if icmpt < ncompartment
        set(h, "facealpha", 0.7);
        set(h, "LineWidth", 1);
    else
        set(h, "LineWidth", 0.03);
    end
    set(h, "facecolor", "interp");
    set(h, "edgealpha", 0.5);

    % Plot compartment numbers
    % (centers approximately computed, without mass matrix)
    center = mean(points, 2);
    center(3) = max(points(3, :)) + 1;
    text(center(1), center(2), center(3), num2str(icmpt), "color", "r", "fontsize", 20);
end
xlabel("x");
ylabel("y");
zlabel("z");
view(3);
grid on;
axis equal;
axis(axis_vec);
title("Finite element mesh, " + subtitle_str);
% title("Finite element mesh", subtitle_str);
