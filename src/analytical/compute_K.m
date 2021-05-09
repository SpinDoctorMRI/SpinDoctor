function K = compute_K(alpha1, n1, alpha2, n2, params)
%COMPUTE_K Compute the integral quantity K.


r = params.r;
D = params.D;
d = params.d;
m = length(D);
R = params.r(m);

dr = r(m) * 1e-12;
r2 = r(1:m) - dr;
r1 = [0 r(1:m-1)] + dr;

K = 0;

if abs(n1 - n2) == 1
    lam1 = alpha1^2;
    lam2 = alpha2^2;
    nu1 = compute_nu(n1, d);
    nu2 = compute_nu(n2, d);
    
    % Computation at point r1
    [V1, dV1] = compute_v(alpha1, n1, r1, params);
    [V2, dV2] = compute_v(alpha2, n2, r1, params);
    I1a = ((lam1 + lam2) * D .* r1.^(d - 1) - (nu1 + nu2 - (d - 1)) * D.^2 .* r1.^(d - 3)) .* V1 .* V2;
    I1b = ((lam2 - lam1) * D .* r1.^d - (nu1 - nu2 - (d - 1)) * D.^2 .* r1.^(d - 2)) .* dV1 .* V2;
    I1c = ((lam1 - lam2) * D .* r1.^d - (nu2 - nu1 - (d - 1)) * D.^2 .* r1.^(d - 2)) .* V1 .* dV2;
    I1d = 2 * D.^2 .* r1.^(d - 1) .* dV1 .* dV2;
    I1 = I1a + I1b + I1c + I1d;
    
    % Computation at point r2
    [V1, dV1] = compute_v(alpha1, n1, r2, params);
    [V2, dV2] = compute_v(alpha2, n2, r2, params);
    I2a = ((lam1 + lam2) * D .* r2.^(d - 1) - (nu1 + nu2 - (d - 1)) * D.^2 .* r2.^(d - 3)) .* V1 .* V2;
    I2b = ((lam2 - lam1) * D .* r2.^d - (nu1 - nu2 - (d - 1)) * D.^2 .* r2.^(d - 2)) .* dV1 .* V2;
    I2c = ((lam1 - lam2) * D .* r2.^d - (nu2 - nu1 - (d - 1)) * D.^2 .* r2.^(d - 2)) .* V1 .* dV2;
    I2d = 2 * D.^2 .* r2.^(d - 1) .* dV1 .* dV2;
    I2 = I2a + I2b + I2c + I2d;
    
    K = sum(I2 - I1) / (lam1 - lam2)^2 / R^(d+1);
end



