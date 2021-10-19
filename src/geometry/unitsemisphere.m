function points = unitsemisphere(npoint)
%UNITSPHERE Create evenly distributed points on unit sphere.
%   This function can be useful to create gradient directions.


% [points, ~, ~] = spheresurface_regularpoints(1, ndirection);
points = fibonacci(2*npoint);
points = points(:, 1:npoint);
end

