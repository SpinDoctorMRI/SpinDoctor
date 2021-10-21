function [j] = J_PGSE(lambda, sequence)
%J_PGSE Compute the quantity J(lambda_n, f) for PGSE sequence
%
%   lambda: double
%   sequence: Sequence class
%
%   j: double


d = sequence.delta;
D = sequence.Delta;

if lambda < 1e-7
    % Use Taylor expansion when lambda is close to 0 
    % to reduce numerical instability
    j = lambda - (D * lambda)^2 / (2 * (D - d/3)) ...
        + ((10*D^3 + 5*D*d^2 - d^3) * lambda^3) / (20 * (3*D - d)) ...
        - ((D^4 + D^2*d^2) * lambda^4) / (8 * (3*D - d)) ...
        + ((21*D^5 + 35*D^3*d^2 + 7*D*d^4 - d^5) * lambda^5) / (840 * (3*D - d));
else
    j = - 1 * ( ...
        + exp(-lambda * (D + d)) ...
        + exp(-lambda * (D - d)) ...
        - 2 * exp(-lambda * d) ...
        - 2 * exp(-lambda * D) ...
        + 2 * (1 - lambda * d)) / ...
        (lambda^2 * d^2 * (D - d/3));
end
end

