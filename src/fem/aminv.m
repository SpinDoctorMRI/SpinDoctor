function [amb, dem] = aminv(ama)
%AMINV Inverse of matrix array.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
%   ama: ama(1:nx, 1:nx, 1:nz)
%   amb: amb(1:nx, 1:nx, 1:nz)

[~, nx, nz] = size(ama);

if nx > 3
    error("Array operation for inverting matrices of dimension more " ...
        + "than 3 (three) is not available, so the performance may be bad.");
end

if nx == 1
    % Matrix determinant
    dem = squeeze(ama);

    % Matrix inverse
    amb = zeros(nx, nx, nz);
    amb(nx, nx, :) = 1 ./ ama(nx, nx, :);
elseif nx == 2

    % Matrix elements
    a = squeeze(ama(1, 1, :)); b = squeeze(ama(1, 2, :));
    c = squeeze(ama(2, 1, :)); d = squeeze(ama(2, 2, :));

    % Matrix determinant
    dem = a .* d - b .* c;

    % Matrix inverse
    amb = zeros(nx, nx, nz);
    amb(1, 1, :) = d ./ dem;
    amb(2, 2, :) = a ./ dem;
    amb(1, 2, :) = -b ./ dem;
    amb(2, 1, :) = -c ./ dem;

else
    % Matrix elements
    a = squeeze(ama(1, 1, :)); b = squeeze(ama(1, 2, :)); c = squeeze(ama(1, 3, :));
    d = squeeze(ama(2, 1, :)); e = squeeze(ama(2, 2, :)); f = squeeze(ama(2, 3, :));
    g = squeeze(ama(3, 1, :)); h = squeeze(ama(3, 2, :)); i = squeeze(ama(3, 3, :));

    % Matrix determinant
    dem = a .* (e .* i - f .* h) - b .* (d .* i - f .* g)+c .* (d .* h - e .* g);

    % Cofactors
    C11 = e .* i - f .* h;
    C12 = -(d .* i - f .* g);
    C13 = d .* h - e .* g;
    C21 = -(b .* i - c .* h);
    C22 = a .* i - c .* g;
    C23 = -(a .* h - b .* g);
    C31 = b .* f - c .* e;
    C32 = -(a .* f - d .* c);
    C33 = a .* e-b .* d;

    % Matrix inverse
    amb = zeros(nx, nx, nz);
    amb(1, 1, :) = C11 ./ dem;
    amb(2, 1, :) = C12 ./ dem;
    amb(3, 1, :) = C13 ./ dem;
    amb(1, 2, :) = C21 ./ dem;
    amb(2, 2, :) = C22 ./ dem;
    amb(3, 2, :) = C23 ./ dem;
    amb(1, 3, :) = C31 ./ dem;
    amb(2, 3, :) = C32 ./ dem;
    amb(3, 3, :) = C33 ./ dem;
end

dem = dem(:)';
