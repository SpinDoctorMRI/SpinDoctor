function points = fibonacci(npoint)
%FIBONACCI Get evenly distributed points on the sphere.
%
%   npoint: int
%
%   points: double(3, npoint)


% Create angles
golden_ratio = (1 + sqrt(5)) / 2;
indices = 0:npoint-1;
theta = 2 * pi * indices / golden_ratio;
phi = acos(1 - 2 * (indices + 0.5) / npoint);

% Create points
points = zeros(3, npoint);
points(1, :) = cos(theta) .* sin(phi);
points(2, :) = sin(theta) .* sin(phi);
points(3, :) = cos(phi);
