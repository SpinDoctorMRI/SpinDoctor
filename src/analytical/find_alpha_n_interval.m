function s = find_alpha_n_interval(xmin, xmax, iteration, n, params)
%FIND_ALPHA_N_INTERVAL Find one root of Fn(alpha) on a small interval interval by dichotomy.
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
%   xmin
%   xmax
%   iteration
%   n
%   params
%
%   s


if iteration < 100
    
    x = (xmin + xmax) / 2;
    s = x;
    
    f = alpha_func(x, n, params);
        
    % if (abs(f)>1e - 8) && (abs(x - xmin) > 1e-13)
    if abs(x - xmin) > 1e-9 % 1e-13
        fmin = alpha_func(xmin, n, params);
        
        if f * fmin < 0
            s = find_alpha_n_interval(xmin, x, iteration + 1, n, params);
        else
            s = find_alpha_n_interval(x, xmax, iteration + 1, n, params);
        end
    end
else
    s = -10;
    warning("Problem in recursive computation: xmin=%e xmax=%e", xmin, xmax);
end
