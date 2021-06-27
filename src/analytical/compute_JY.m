function [J, Y, dJ, dY] = compute_JY(z, n, d)
%COMPUTE_JY Compute the Bessel functions J and Y.
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
%   z
%   n
%   d
%
%   J
%   Y
%   dJ
%   dY


switch d
    case 2
        J = besselj(n, z);
        Y = bessely(n, z);
        dJ = (besselj(n - 1, z) - besselj(n + 1, z)) / 2;
        dY = (bessely(n - 1, z) - bessely(n + 1, z)) / 2;
    case 3
        J = besseljs(n, z);
        Y = besselys(n, z);
        
        ind = find(z > 0);
        zz = z(ind);
        dJ(ind) = (besseljs(n - 1, zz) - besseljs(n + 1, zz) - besseljs(n, zz) ./ zz ) / 2;
        dY(ind) = (besselys(n - 1, zz) - besselys(n + 1, zz) - besselys(n, zz) ./ zz ) / 2;
        
        ind = find(z == 0);
        dJ(ind) = (n == 1) / 3;
        dY(ind) = Inf;               % check the sign!
end
