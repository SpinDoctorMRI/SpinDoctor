function plot_mesh(femesh)


% Determine limits of domain
pmin = min([femesh.points{:}], [], 2);
pmax = max([femesh.points{:}], [], 2);
axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];
facets = [femesh.facets{1, :}];
points = femesh.points{1};
h = trisurf(facets', points(1, :), points(2, :), points(3, :));
set(h, "LineWidth", 0.01);
set(h, "facecolor", "interp");
set(h, "edgealpha", 0.15);
xlabel("x");
ylabel("y");
zlabel("z");
view(3);
grid on;
axis equal;
axis(axis_vec);
title("Finite element mesh");