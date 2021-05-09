function [alpha, n] = compute_eigenvalues(alpha_max, params, dalpha)
%COMPUTE_EIGENVALUES Compute radial and angular eigenvalues.
%   Find all radial eigenvalues smaller than alpha_max^2, with radial and angular indices.
%
% dalpha – smallest distance between two alphas
% alpha  - sqrt lambda (radial eigenvalue)
% n      - angular parameter: n = 0 (1D), n^2 = nu (2D), or n(n+1) = nu (3D)


alpha = [];
n = [];

alpha0 = 0;
j = 0;
while ~isempty(alpha0) && alpha0(1) < alpha_max
    fprintf('Computing zeros of F(alpha) for n = %d\n', j);
    alpha0 = find_alpha_n(alpha_max, j, params, dalpha);
    alpha = [alpha alpha0];
    n = [n j*ones(size(alpha0))];
    j = j + 1;
end

[alpha, inds] = sort(alpha);
n = n(inds);


% Number of layers
L = length(params.D);

%if params.W(length(params.D)) < 1e-12 && ((params.D(1) > 1e-12) || ((params.D(1) < 1e-12) && (params.W(end) < 1e-12)) )
if params.W(L) < 1e-12 && (length(params.W) == L || params.W(end) < 1e-12)
    % Add ground eigenmode explicitly
    alpha = [0 alpha];
    n = [0 n];
end

inds_keep = alpha < alpha_max;
alpha = alpha(inds_keep);
n = n(inds_keep);
