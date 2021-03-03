function amb = amtam(amx, ama)
%AMTAM Matrix array transpose times matrix array.
%   Copyright (c) 2013,  Talal Rahman,  Jan Valdman
%
%   ama: ama(1:nx, 1:ny, 1:nz)
%   amx: amx(1:nx, 1:nk, 1:nz)
%   amb: amb(1:nk, 1:ny, 1:nz)

[~, ny, nz] = size(ama);
[~, nk, ~] = size(amx);

amb = zeros(nk, ny, nz);
for row = 1:nk
    amb(row, :, :) = avtam(amx(:, row, :), ama);
end
