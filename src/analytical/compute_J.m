function J = compute_J(alpha, params)
%COMPUTE_J Compute the integral quantity J.
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
%   alpha
%   params
%
%   J


r = params.r;
D = params.D;
W = params.W;
d = params.d;

m = length(D);
R = r(m);

lam = alpha^2;

if alpha > 1e-6
    V = compute_v(alpha, 0, R, params);
    J = W(m) * V / lam / R;
else
    J = 1 / d;
end
