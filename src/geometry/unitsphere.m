function points = unitsphere(npoint)
%UNITSPHERE Create evenly distributed points on unit sphere.
%   This function can be useful to create gradient directions.


% [points, ~, ~] = spheresurface_regularpoints(1, ndirection);
points = fibonacci(npoint);

% remove negative zeros
points(points == 0) = +0;

end

