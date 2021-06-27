function avb = svamt(svx, ama)
%SVAMT Scalar vector times array of matrices transposed.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   svx: [1 x nx]
%   ama: [ny x nx x nz]
%
%   avb: [1 x ny x nz]


[ny, ~, nz] = size(ama);

avx = svx;
avx = avx(ones(ny, 1), :, ones(nz, 1));

avb = ama .* avx;
avb = sum(avb, 2);
avb = reshape(avb, 1, ny, nz);
