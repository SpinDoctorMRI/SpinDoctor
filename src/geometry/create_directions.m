function directions = create_directions(ndirection, flat_dirs, remove_opposite)
%CREATE_DIRECTIONS Find directions distributed on the sphere or in the plane.
%
%   ndirection:
%       Number of points on the unit sphere to create.
%   flat_dirs: logical
%         Indicator for points on unit circle instead of unit sphere.
%   remove_opposite: logical
%       Indicator for removing redundant (opposite) directions.
%
%   directions: struct with fields
%       points: double(3, npoints)
%           Gradient directions.
%       indices: int(1, ndir_unique)
%           Indices of the directions to keep.
%       opposite: cell(1, ndir_unique)
%           Indices of opposite directions.
%
% where ndir_unique is the number of unique directions, and ndirection is the total
% number of directions, including opposite directions.
%
% points(:, i) are the coordinates of direction i.
% indices(j) is the index of the j'th direction to be kept.
% opposite{j} is equal to:
%     [] if the direction indices(j) has no opposite direction
%     i  if the direction indices(j) has an opposite direction i, and in
%         that case points(:, i) = - points(:, indices(j))
%
% Example: [1 2 3 4 5 6 7 8] are the total number of directions, with direction
% 7 opposite of direction 5, and direction 3 opposite of direction 2. Then:
% ndirection = 8
% ndir_unique = 6
% indices = [1 2 4 5 6 8]
% opposite = {[] 3 [] 7 [] []}


% Tolerance for norm equality between opposite directions
tolerance = 1e-10;

% Create points
if flat_dirs
    % Angles
    thetavec = linspace(0, 2 * pi, ndirection + 1)';
    % thetavec = linspace(0, pi, ndirection + 1)';

    % Create points on the unit circle in the x-y plane
    points = zeros(3, ndirection);
    points(1, :) = cos(thetavec(1:end-1));
    points(2, :) = sin(thetavec(1:end-1));
else
    % Create points on the unit sphere
    points = fibonacci(ndirection);
    % [points, ~, ~] = spheresurface_regularpoints(1, ndirection);
    ndirection = size(points, 2);
end

% Initialize output arguments
indices = 1:ndirection;
opposite = cell(1, ndirection);

% Identify opposite directions (assuming no direction duplicates)
if remove_opposite
    for ikeep = 1:ndirection-1
        for iremove = ikeep+1:ndirection
            % Check if the two directions are opposite (up to a tolerance). We
            % assume the points are of unit norm, so their scalar product has a
            % lower bound of -1
            if points(:, ikeep)' * points(:, iremove) + 1 < tolerance
                % Mark indices to be removed with 0
                indices(iremove) = 0;

                % Associate kept with removed index
                opposite{ikeep} = iremove;
            end
        end
    end

    % Remove indices
    indices(indices == 0) = [];
end

% % Create output structure
% directions.points = points;
% directions.indices = indices;
% directions.opposite = opposite;
% return;

% ipoints = 1:8;
% indices = [1 2 4 5 6 8];
% opposite = {[] 3 [] 7 [] []};

% Reorder points so that indices are monotonically increasing integer values.
% This is only needed for parfor loops to be able to loop through indices
indices_reorder = 1:length(indices);
ipoint_reorder = [indices, opposite{:}];
[~, ipoint_reorder_inverse] = sort(ipoint_reorder);
opposite_reorder = cellfun(@(i) ipoint_reorder_inverse(i), opposite, "UniformOutput", false);
points_reorder = points(:, ipoint_reorder);

% Create output structure
directions.points = points_reorder;
directions.indices = indices_reorder;
directions.opposite = opposite_reorder;
