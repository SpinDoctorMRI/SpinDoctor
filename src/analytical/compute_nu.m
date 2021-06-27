function nu = compute_nu(n, d)
%COMPUTE_NU Compute the angular eigenvalue nu as a function of the angular index n.
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
%   n
%   d
%
%   nu


switch d
  case 2
    nu = n^2;
  case 3
    nu = n * (n + 1); 
end





