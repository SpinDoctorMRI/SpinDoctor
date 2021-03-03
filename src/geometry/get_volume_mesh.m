function [total_volume, volumes, centers] = get_volume_mesh(points, elements)
%GET_VOLUME_MESH Compute volumes and centers of finite elements.
%
%   points: double(3, npoints)
%   elements: int(4, nelement)
%
%   total_volume: double
%   volumes: double(1, nelement)
%   centers: double(3, nelement)


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
