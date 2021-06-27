function amb = smamt(smx, ama)
%SMAMT Scalar matrix times array of matrices transposed.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   smx: [nk x nx]
%   ama: [ny x nx x nz]
%
%   amb: [nk x ny x nz]


[ny, ~, nz] = size(ama);
[nk, ~] = size(smx);

amb = zeros(nk, ny, nz);
for row = 1:nk
    amb(row, :, :) = svamt(smx(row, :), ama);
end
