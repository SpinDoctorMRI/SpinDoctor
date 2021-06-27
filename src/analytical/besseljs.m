function js = besseljs(nu, x)
%BESSELJS Sperical Bessel function jnu(x)
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
%   nu
%   x
%
%   js


lnu = size(nu, 2);

xm = repmat(x, 1, lnu);
js = sqrt(pi ./(2 * xm)) .* besselj(nu + 0.5, x);
