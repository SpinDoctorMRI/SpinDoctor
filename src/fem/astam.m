function amb = astam(asx, ama)
% Copyright (c) 2015, Talal Rahman, Jan Valdman
%
%   ama: [nx x ny x nz]
%   asx: [1 x nz]
%
%   amb: [nx x ny x nz]

[nx, ny, nz] = size(ama);

asx = reshape(asx, 1, 1, nz);
asx = asx(ones(nx, 1), ones(ny, 1), :);

amb = ama .* asx;
