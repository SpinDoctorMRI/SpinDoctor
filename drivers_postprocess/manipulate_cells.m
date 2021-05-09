cells = create_cells(setup);


%%
cells = add_cyl( ...
cells, [
-25.5
21.4
], [
2
]);

%%
cells.centers(:, [17]) = [
0.1443
2.5826
];

%%
cells.radii([17]) = [
2.4
];

%%
cells = touch(setup, cells, 17);

%%
cells = rem_cyl(cells, []);



%% Center
mass_center = cells.centers * cells.radii'.^2 / sum(cells.radii.^2);
cells.centers = cells.centers - mass_center;

%% Rotate
theta = -pi / 40;
R = [
    cos(theta) sin(theta)
    -sin(theta) cos(theta)
];
cells.centers = R * cells.centers;

%% Plot configuration
close all
% plot_cells(cells, setup)
plot_cells_2(cells, setup)

%%
function cells = add_cyl(cells, c, r)
    cells.centers = [cells.centers c];
    cells.radii = [cells.radii r];
end

function cells = rem_cyl(cells, icell)
    cells.centers(:, icell) = [];
    cells.radii(icell) = [];
end

function cells = touch(setup, cells, icell)
    rmean = (setup.geometry.rmin + setup.geometry.rmax) / 2;
    ncell = length(cells.radii);
    
    dist = vecnorm(cells.centers - cells.centers(:, icell)) - cells.radii;
    dist = min(dist([1:icell-1 icell+1:end]));
    rmin1 = dist - setup.geometry.dmax * rmean;
    rmax1 = dist - setup.geometry.dmin * rmean;
    r_lower = max(setup.geometry.rmin, rmin1);
    r_upper = min(setup.geometry.rmax, rmax1);
    if r_lower <= r_upper
        radius = (r_lower + r_upper) / 2;
        % radius = r_upper;
        cells.radii(icell) = radius;
    else
        error("Uncompatible");
    end
end
