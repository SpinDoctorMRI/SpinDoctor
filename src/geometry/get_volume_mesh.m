function [total_volume, volumes, centers] = get_volume_mesh(points, elements)
%GET_VOLUME_MESH Compute volumes and centers of finite elements.
%
%   points: [3 x npoint]
%   elements: [4 x nelement]
%
%   total_volume: [1 x 1]
%   volumes: [1 x nelement]
%   centers: [3 x nelement]


% Sizes
nelement = size(elements, 2);

% Elements
x = reshape(points(:, elements), 3, 4, nelement);

% Element centers
centers = squeeze(mean(x, 2));

% Element volumes
areavectors = cross(x(:, 2, :) - x(:, 4, :), x(:, 3, :) - x(:, 4, :));
volumes = 1 / 6 * shiftdim(abs(dot(x(:, 1, :) - x(:, 4, :), areavectors)), 1);

% Total volume
total_volume = sum(volumes);
