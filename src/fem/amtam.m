function amb = amtam(amx, ama)
%AMTAM Matrix array transpose times matrix array.
%   Copyright (c) 2013,  Talal Rahman,  Jan Valdman
%
%   amx: [nx x nk x nz]
%   ama: [nx x ny x nz]
%
%   amb: [nk x ny x nz]

[~, ny, nz] = size(ama);
[~, nk, ~] = size(amx);

amb = zeros(nk, ny, nz);
for row = 1:nk
    amb(row, :, :) = avtam(amx(:, row, :), ama);
end
