function plot_cells_2d(cells, setup)
%PLOT_CELLS_2D Plot cells of canonical configuration in the plane.
%
%   cells: struct
%   setup: struct

centers = cells.centers;
radii = cells.radii;

ncell = size(centers, 2);

figure;
hold on

% Cylinders
for icell = 1:ncell
    linewidth = 0.5;
    if 200 < icell
        linewidth = 1.5;
    end
    [X, Y] = cylinder;
    X = radii(icell) * X + centers(1, icell);
    Y = radii(icell) * Y + centers(2, icell);
    plot(X(1, :), Y(1, :), "linewidth", linewidth);
    p = centers(:, icell);
    text(p(1), p(2), num2str(icell), "color", "r", "fontsize", 8);
end

xlabel("x");
ylabel("y");
grid on
axis equal
axis tight
title("Cells of canonical configuration");
