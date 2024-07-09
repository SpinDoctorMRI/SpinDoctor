function points = unitsemisphere(npoint, normal)
%UNITSPHERE Create evenly distributed points on unit sphere.
%   This function can be useful to create gradient directions.


% [points, ~, ~] = spheresurface_regularpoints(1, ndirection);
points = fibonacci(2*npoint);
points = points(:, 1:npoint);

if nargin == nargin(@unitsemisphere) && (normal(1) ~= 0 || normal(2) ~= 0)
    % Create rotation matrix to transform normal to e_z
    normal = normal / norm(normal);
    c = normal + [0; 0; 1];
    R = 2  * (c * c') / (c' * c) - eye(3);
    points = R' * points;
end

% remove negative zeros
points(points == 0) = +0;

end

