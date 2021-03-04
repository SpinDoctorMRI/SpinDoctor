function results = compute_blochtorrey_eig(q, lap_eig, directions)
%COMPUTE_BLOCHTORREY_EIG Compute the eigenvectors and eigenvalues of the Bloch-Torrey operator.
%
%   q: double
%       q-value
%      lap_eig: struct with fields
%       values: double(neig, 1)
%       moments: double(neig, neig, 3)
%     directions: struct
%       Gradient directions.
%
%   results: struct with fields
%        Vsort: cell(1, ndirection)
%       Dsort: cell(1, ndirection)
%        invVsortC1: cell(1, ndirection)
%        invVsort: cell(1, ndirection)
%         totaltime: double


% Measure computational time
starttime = tic;

% Extract HARDI points
dir_points = directions.points;
dir_inds = directions.indices;
% dir_opposite = directions.opposite;

% Sizes
ndirection = size(dir_points, 2);

% Initialize output arguments
Vsort = cell(1, ndirection);
Dsort = cell(1, ndirection);
invVsort = cell(1, ndirection);
invVsortC1 = cell(1, ndirection);

% Perform eigendecomposition for each compartment and direction
L_mat = diag(lap_eig.values);
for idir = dir_inds
    gdir = dir_points(:, idir);
    W_mat = sum(lap_eig.moments .* shiftdim(gdir, -2), 3);
    K_mat = L_mat + 1i * W_mat * q;
    [V, D] = eig(K_mat);
    [~, inds_sort] = sort(real(diag(D)));
    V = V(:, inds_sort);
    D = diag(D(inds_sort, inds_sort));
    invV = inv(V);
    c1 = zeros(size(V, 1), 1);
    c1(1) = 1;
    invVC1 = V \ c1;
    Vsort{idir} = V;
    Dsort{idir} = D;
    invVsort{idir} = invV;
    invVsortC1{idir} = invVC1;
end

% Create output structure
results.Vsort = Vsort;
results.Dsort = Dsort;
results.invVsort = invVsort;
results.invVsortC1 = invVsortC1;
results.totaltime = toc(starttime);

% Display time of function evaluation
toc(starttime)
