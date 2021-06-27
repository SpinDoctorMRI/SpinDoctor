function avb = avtam(avx, ama)
%AVTAM Vector array transpose times matrix array.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   avx: [nx x nz]
%   ama: [nx x ny x nz]
%
%   avb: [1 x ny xnz]

[~, ny, ~] = size(ama);

avx = avx(:, ones(ny, 1), :);

avb = ama .* avx;
avb = sum(avb, 1);
