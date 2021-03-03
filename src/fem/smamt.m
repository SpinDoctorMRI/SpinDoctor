function amb = smamt(smx, ama)
%SMAMT Scalar matrix times array of matrices transposed.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   ama: ama(1:ny, 1:nx, 1:nz)
%   smx: smx(1:nk, 1:nx)
%   amb: amb(1:nk, 1:ny, 1:nz)

[ny, ~, nz] = size(ama);
[nk, ~] = size(smx);

amb = zeros(nk, ny, nz);
for row = 1:nk
    amb(row, :, :) = svamt(smx(row, :), ama);
end
