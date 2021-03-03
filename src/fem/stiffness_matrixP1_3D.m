function [A, volumes] = stiffness_matrixP1_3D(elements, coordinates, coeffs)
%STIFFNESS_MATRIXP1_3D Assemble 3D stiffness matrix using P1 elements.
%   Copyright (c) 2016, Jan Valdman
%
% coeffs can be either P0 (elementwise constant) or P1 (elementwise nodal)
% function represented by a collumn vector with size(elements, 1) or
% size(coordinates, 1) entries if coeffs is not provided then coeffs = 1 is
% assumed globally)


% Number of elements
NE = size(elements, 1);

% Problem dimension
DIM = size(coordinates, 2);

% Number of local basic functions
NLB = 4;

% Particular part for a given element in a given dimension
coord = zeros(DIM, NLB, NE);
for d = 1:DIM
    for i = 1:NLB
        coord(d, i, :) = coordinates(elements(:, i), d);
    end
end
IP = [1/4 1/4 1/4]';
[dphi, jac] = phider(coord, IP, "P1");

volumes = abs(squeeze(jac)) / factorial(DIM); % det(J)/3

dphi = squeeze(dphi);

if nargin < 3
    Z = astam(volumes', amtam(dphi, dphi));
elseif isvector(coeffs)
    if numel(coeffs) == size(coordinates, 1)
        % P1->P0 averaging
        % coeffs = evaluate_average_point(elements, coeffs);
        error("Not implemented");
    end
    Z = astam((volumes .* coeffs)', amtam(dphi, dphi));
else
    % Coeffs is either a scalar tensor (3 x 3) or an (3 x 3 x nelement) or (3 x
    % 3 x nnode)
    if size(coeffs, 3) == 1
        % Apply diffusion tensor to dphi, then dphi' to the results. Weigh on
        % each element by volumes
        Z = astam(volumes', amtam(dphi, smamt(coeffs, amt(dphi))));
    elseif size(coeffs, 3) == size(coordinates, 1)
        % size(coeffs) = 3 x 3 x nnode
        % each tensor in coeffs is assumed symmetric
        % average points tensors in each element
        error("Not implemented");
    else
        % size(coeffs) = 3 x 3 x nelement
        % each tensor in coeffs is assumed symmetric
        Z = astam(volumes', amtam(dphi, amtam(amt(coeffs), dphi)));
    end
end

Y = reshape(repmat(elements, 1, NLB)', NLB, NLB, NE);

X = permute(Y, [2 1 3]);
A = sparse(X(:), Y(:), Z(:));

% MODFICATION: Assure symmetry
sym = @(x) (x + x') / 2;
A = sym(A);
