function points = unitsemicircle(npoint, normal)
%UNITSEMICIRCLE Create evenly distributed points on unit semicircle.
%   This function can be useful to create gradient directions.
%
%   UNITSEMICIRCLE(NPOINT) creates NPOINT evenly distributed points in the x-y
%   halfplane (with y >= 0).
%
%   UNITSEMICIRCLE(NPOINT, NORMAL) distributes the point in the plane orthogonal to
%   NORMAL.
%
%   npoint: [1 x 1]
%   normal (optional): [3 x 1]
%
%   points: [3 x npoint]


% Angles
angles = linspace(0, pi, ndirection + 1)';

% Create points on the unit circle in the x-y plane
points = zeros(3, npoint);
points(1, :) = cos(angles(1:end-1));
points(2, :) = sin(angles(1:end-1));

if nargin == nargin(@unitsemicircle) && (normal(1) ~= 0 || normal(2) ~= 0)
    % Create rotation matrix to transform normal to e_z
    normal = normal / norm(normal);
    c = normal + [0; 0; 1];
    R = 2  * (c * c') / (c' * c) - eye(3);
    points = R' * points;
end

end

