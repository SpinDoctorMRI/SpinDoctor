function points = unitcircle(npoint, normal)
%UNITCIRCLE Create evenly distributed points on unit circle.
%   This function can be useful to create gradient directions.
%
%   UNITCIRCLE(NPOINT) creates NPOINT evenly distributed points in the x-y
%   plane.
%
%   UNITCIRCLE(NPOINT, NORMAL) distributes the point in the plane orthogonal to
%   NORMAL.
%
%   npoint: [1 x 1]
%   normal (optional): [3 x 1]
%
%   points: [3 x npoint]


% Angles
angles = linspace(0, 2 * pi, npoint + 1)';

% Create points on the unit circle in the x-y plane
points = zeros(3, npoint);
points(1, :) = cos(angles(1:end-1));
points(2, :) = sin(angles(1:end-1));

if nargin == nargin(@unitcircle) && (normal(1) ~= 0 || normal(2) ~= 0)
    % Create rotation matrix to transform normal to e_z
    normal = normal / norm(normal);
    c = normal + [0; 0; 1];
    R = 2  * (c * c') / (c' * c) - eye(3);
    points = R' * points;
end

% remove negative zeros
points(points == 0) = +0;

end

