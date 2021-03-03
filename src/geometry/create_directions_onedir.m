function directions = create_directions_onedir(direction)
%CREATE_DIRECTIONS_ONEDIR Create one direction with no opposites.
%
%   direction: double(3, 1)
%       Gradient direction.
%
%   directions: struct with fields
%       points: double(3, npoints)
%           Gradient directions.
%       indices: int(1, ndir_unique)
%           Indices of the directions to keep.
%       opposite: cell(1, ndir_unique)
%           Indices of opposite directions.


% Create output structure
directions.points = direction / norm(direction);
directions.indices = 1;
directions.opposite = {[]};
