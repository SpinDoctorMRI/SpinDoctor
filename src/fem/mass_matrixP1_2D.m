function M = mass_matrixP1_2D(elements, areas, coeffs)
%MASS_MATRIXP1_2D Assemble 2D mass matrix using P1 elements.
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
% coeffs can be only P0 (elementwise constant) function
% represented by a collumn vector with size(elements, 1) entries
% if coeffs is not provided then coeffs = 1 is assumed globally
% Note: P1 coeffs needs a higher integration rule (not implemented yet)


X = kron(ones(1, 3), elements);
Y = kron(elements, ones(1, 3));

if nargin < 3
    Z = kron(areas, reshape((ones(3) + eye(3)) / 12, 1, 9));
else
    if numel(coeffs) == size(elements, 1)
        % P0 coefficients
        Z = kron(areas .* coeffs, reshape((ones(3)+eye(3)) / 12, 1, 9));
    else
        % P1 coefficients
        M1 = [6 2 2; 2 2 1; 2 1 2] / 60;
        M2 = M1([3, 1, 2], [3, 1, 2]);
        M3 = M2([3, 1, 2], [3, 1, 2]);

        Z = kron(areas .* coeffs(elements(:, 1)), reshape(M1, 1, 9)) ...
            + kron(areas .* coeffs(elements(:, 2)), reshape(M2, 1, 9)) ...
            + kron(areas .* coeffs(elements(:, 3)), reshape(M3, 1, 9));
    end
end

M = sparse(X, Y, Z); % M(X(k), Y(k)) = Z(k)
