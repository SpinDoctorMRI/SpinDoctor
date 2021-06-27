function I = compute_I(alpha1, alpha2, n, params)
%COMPUTE_I Compute the integral quantity I.
%
% This function is based on the following articles and corresponding code:
%   [1] D. S. Grebenkov, NMR Survey of Reflected Brownian Motion,
%       Rev. Mod.Phys. 79, 1077 (2007)
%   [2] D. S. Grebenkov, Pulsed-gradient spin-echo monitoring of restricted 
%       diffusion inmultilayered structures,
%       J. Magn. Reson. 205, 181-195 (2010).
%
%   alpha1
%   alpha2
%   n
%   params
%
%   I


r = params.r;
D = params.D;
d = params.d;

m = length(D);
R = params.r(m);

dr = r(m) * 1e-12;
r2 = r(1:m) - dr;
r1 = [0 r(1:m-1)] + dr;

% Check if alpha1 is not equal to alpha2
if abs(alpha1 - alpha2) > 1e-8
    lam1 = alpha1^2;
    lam2 = alpha2^2;
    
    % Computation at point r1
    [V1, dV1] = compute_v(alpha1, n, r1, params);
    [V2, dV2] = compute_v(alpha2, n, r1, params);
    I1 = D .* r1.^(d - 1) .* (dV1 .* V2 - dV2 .* V1);
    
    % Computation at point r2
    [V1, dV1] = compute_v(alpha1, n, r2, params);
    [V2, dV2] = compute_v(alpha2, n, r2, params);
    I2 = D .* r2.^(d - 1) .* (dV1 .* V2 - dV2 .* V1);
    
    I = (I2 - I1) / (lam1 - lam2) / R^d;
else
    alpha = alpha1;
    lam = alpha^2;
    nu = compute_nu(n, d);
    
    if alpha > 1e-6
        [V, dV] = compute_v(alpha, n, r1, params);
        I1 = D .* r1.^d .* dV.^2 ...
            + (lam * r1.^d - D .* r1.^(d - 2) * nu) .* V.^2 ...
            + (d - 2) * D .* r1.^(d - 1) .* dV .* V;
        
        [V, dV] = compute_v(alpha, n, r2, params);
        I2 = D .* r2.^d .* dV.^2 ...
            + (lam * r2.^d - D .* r2.^(d - 2) * nu) .* V.^2 ...
            + (d - 2) * D .* r2.^(d - 1) .* dV .* V;
        
        I = (I2 - I1) / (2 * lam * R^d);
    else
        I = (r2.^d - r1.^d) / R^d / d;
    end
end



