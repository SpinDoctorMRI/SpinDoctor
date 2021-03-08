function cells = create_cells(setup)
%CREATE_CELLS Create random cell centers and radii for cylinders or spheres.
%
%   params_cell: struct
%
%   cells: struct with fields
%       centers: [d x ncell], where d is 2 (cylinders) or 3 (spheres)
%       radii: [1 x ncell]


% Extract cell parameters
cell_shape = setup.geometry.cell_shape;
ncell = setup.geometry.ncell;
rmin = setup.geometry.rmin;
rmax = setup.geometry.rmax;
dmin = setup.geometry.dmin;
dmax = setup.geometry.dmax;

% Spatial dimension
if cell_shape == "sphere"
    d = 3;
elseif cell_shape == "cylinder"
    d = 2;
end

% Maximum number of random cell centers to try out
npoint_max = 100000 * ncell;

% Mean allowed radius
rmean = (rmin + rmax) / 2;

% Cell centers (first is zero)
centers = zeros(d, ncell);

% Cell radii (first is rmean)
radii = zeros(1, ncell);
radii(1) = rmean;

% Initiate iterators
icell = 1;
npoint = 1;

% Generate random cell centers until ncell cells are created
while icell < ncell && npoint < npoint_max

    % Generate a random point using a uniform distribution with zero mean and
    % variance proportional to rmean
    if cell_shape == "sphere"
        point = (rand(3, 1) - 0.5) * rmean * max(10, ncell^(1 / 3));
    else
        point = (rand(2, 1) - 0.5) * rmean * max(10, sqrt(ncell)) * 40;
    end

    % Distance from point to the (already determined) cells
    dist = vecnorm(centers - point) - radii;
    dist = dist(1:icell);

    % Distance from point to nearest cell
    dist = min(dist);

    % Maximum allowed distance to nearest cell gives lowest allowed radius
    rmin1 = dist - dmax * rmean;

    % Minimum allowed distance to nearest cell gives highest allowed radius
    rmax1 = dist - dmin * rmean;

    % Radius has to be higher than both of the lowest allowed radii
    r_lower = max(rmin, rmin1);

    % Radius has to be lower than both of the highest allowed radii
    r_upper = min(rmax, rmax1);

    % Check if any radius is allowed
    if r_lower <= r_upper
        % Take midpoint between the bounds of all allowed radii
        radius = (r_lower + r_upper) / 2;

        % Update cells
        icell = icell + 1;
        centers(:, icell) = point;
        radii(icell) = radius;
    end

    % Update count for number of random centers generated
    npoint = npoint + 1;
end

% Check that all cells were created
if icell < ncell
    error("Did not find enough cell centers. Consider using larger separation tolerances.");
end

% Center cell collection
pmin = min(centers - radii, [], 2);
pmax = max(centers + radii, [], 2);
pmean = (pmin + pmax) / 2;
centers = centers - pmean;

% Create output structure
cells.centers = centers;
cells.radii = radii;
