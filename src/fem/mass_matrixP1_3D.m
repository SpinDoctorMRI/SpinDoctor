function M = mass_matrixP1_3D(elements, volumes, coeffs)
%MASS_MATRIXP1_2D Assemble 3D mass matrix using P1 elements.
%   Copyright (c) 2016, Jan Valdman
%
% coeffs can be only P0 (elementwise constant) function
% represented by a collumn vector with size(elements, 1) entries
% if coeffs is not provided then coeffs = 1 is assumed globally
% Note: P1 coeffs needs a higher integration rule (not implemented yet)

X = kron(ones(1, 4), elements);
Y = kron(elements, ones(1, 4));

if nargin < 3
    Z = kron(volumes, reshape((ones(4) + eye(4)) / 20, 1, 16));
else
    if numel(coeffs) == size(elements, 1)
        % P0 coefficients
        Z = kron(volumes .* coeffs, reshape((ones(4) + eye(4)) / 20, 1, 16));

    else
        % P1 coefficients
        M1 = [6 2 2 2; 2 2 1 1; 2 1 2 1; 2 1 1 2]/120;
        M2 = M1([4, 1, 2, 3], [4, 1, 2, 3]);
        M3 = M2([4, 1, 2, 3], [4, 1, 2, 3]);
        M4 = M3([4, 1, 2, 3], [4, 1, 2, 3]);

        Z = kron(volumes .* coeffs(elements(:, 1)), reshape(M1, 1, 16)) ...
            +kron(volumes .* coeffs(elements(:, 2)), reshape(M2, 1, 16)) ...
            +kron(volumes .* coeffs(elements(:, 3)), reshape(M3, 1, 16)) ...
            +kron(volumes .* coeffs(elements(:, 4)), reshape(M4, 1, 16));
    end
end

M = sparse(X, Y, Z);

% MODFICATION: Enforce symmetry
sym = @(x) (x + x') / 2;
M = sym(M);
