function amb = amsm(ama, smx)
%AMSM Matrix array times scalar matrix.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   ama: [nx x ny x nz]
%   smx: [ny x nk]
%
%   amb: [nx x nk x nz]

[nx, ~, nz] = size(ama);
[~, nk] = size(smx);

amb = zeros(nx, nk, nz);
for col = 1:nk
    amb(:, col, :) = amsv(ama, smx(:, col));
end
