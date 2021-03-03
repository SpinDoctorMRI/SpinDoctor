function amb = amsm(ama, smx)
%AMSM Matrix array times scalar matrix.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   ama: ama(1:nx, 1:ny, 1:nz)
%   smx: smx(1:ny, 1:nk)
%   amb: amb(1:nx, 1:nk, 1:nz)

[nx, ~, nz] = size(ama);
[~, nk] = size(smx);

amb = zeros(nx, nk, nz);
for col = 1:nk
    amb(:, col, :) = amsv(ama, smx(:, col));
end
