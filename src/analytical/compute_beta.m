function beta = compute_beta(alpha, n, params)
%COMPUTE_BETA Compute normalization coefficients beta.
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
%   alpha
%   n
%   params
%
%   beta


I = compute_I(alpha, alpha, n, params);
beta = 1 / sqrt(2 * sum(I));
