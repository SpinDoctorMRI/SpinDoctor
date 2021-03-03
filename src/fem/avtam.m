function avb = avtam(avx, ama)
%AVTAM Vector array transpose times matrix array.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   ama: ama(1:nx, 1:ny, 1:nz)
%   avx: avx(1:nx, 1:nz)
%   avb: avb(1, 1:ny, 1:nz)

[~, ny, ~] = size(ama);

avx = avx(:, ones(ny, 1), :);

avb = ama .* avx;
avb = sum(avb, 1);
