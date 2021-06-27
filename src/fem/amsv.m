function avb = amsv(ama, svx)
%AMSV Matrix array times scalar vector.
%   Copyright (c) 2013,  Talal Rahman,  Jan Valdman
%
%   ama: [nx x ny x nz]
%   svx: [ny x 1]
%
%   avb: [nx x 1 x nz]

[nx, ~, nz] = size(ama);

avx = svx(:).';
avx = avx(ones(nx, 1), :, ones(nz, 1));

avb = ama .* avx;
avb = sum(avb, 2);
