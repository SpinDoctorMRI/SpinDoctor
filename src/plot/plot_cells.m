function plot_cells(cells, setup)
%PLOT_CELLS Plot cells of canonical configuration.
%
%   cells: struct
%   setup: struct

centers = cells.centers;
radii = cells.radii;

[d, ncell] = size(centers);

figure;
hold on;
for icell = 1:ncell
    if d == 2
        % Cylinders
        [X, Y, Z] = cylinder;
        X = radii(icell) * X + centers(1, icell);
        Y = radii(icell) * Y + centers(2, icell);
        Z = setup.geometry.height * (Z - 0.5);
        surf(X, Y, Z);%, "edgealpha", 0);
        p = [centers(:, icell); setup.geometry.height / 2 + 1];
    else
        % Spheres
        [X, Y, Z] = sphere;
        X = radii(icell) * X + centers(1, icell);
        Y = radii(icell) * Y + centers(2, icell);
        Z = radii(icell) * Z + centers(3, icell);
        surf(X, Y, Z);%, "facealpha", 0.6, "edgealpha", 0.5);
        p = centers(:, icell) + [0; 0; 1] * (radii(icell) + 1);
    end
    text(p(1), p(2), p(3), num2str(icell), "color", "r", "fontsize", 20);
end

xlabel("x");
ylabel("y");
zlabel("z");
grid on
axis equal
view(3);
title("Cells of canonical configuration");
