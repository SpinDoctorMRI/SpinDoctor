function beta = compute_beta(alpha, n, params)
%COMPUTE_BETA Compute normalization coefficients beta.


I = compute_I(alpha, alpha, n, params);
beta = 1 / sqrt(2 * sum(I));
