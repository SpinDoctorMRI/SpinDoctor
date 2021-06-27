function bc = compute_bc(alpha, n, params)
%COMPUTE_BC Compute the coefficients b and c.
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
%   bc


r = params.r;
D = params.D;
W = params.W;
d = params.d;
m = length(D);

bc = zeros(m, 2);
bc(:, 1) = 1;

if alpha == 0
    return;
end

rho = r(1:m-1) ./ sqrt(D(1:m-1));
rhot = r(1:m-1) ./ sqrt(D(2:m));

for i = 1:m-1
    [J, Y, dJ, dY] = compute_JY(alpha * rho(i), n, d);
    [J1, Y1, dJ1, dY1] = compute_JY(alpha * rhot(i), n, d);

    A11 = -sqrt(D(i) / D(i+1)) * dJ * Y1 + J * dY1 + sqrt(D(i)) / W(i) * alpha * dJ * dY1;
    A12 = -sqrt(D(i) / D(i+1)) * dY * Y1 + Y * dY1 + sqrt(D(i)) / W(i) * alpha * dY * dY1;
    A21 = sqrt(D(i) / D(i+1)) * dJ * J1 - J * dJ1 - sqrt(D(i)) / W(i) * alpha * dJ * dJ1;
    A22 = sqrt(D(i) / D(i+1)) * dY * J1 - Y * dJ1 - sqrt(D(i)) / W(i) * alpha * dY * dJ1;

    z = alpha * rhot(i);
    switch d
        case 2
            Q = z * pi / 2;
        case 3
            Q = z^2;
    end
    A = Q * [A11 A12; A21 A22];
    bc(i + 1, :) = A * bc(i, :)';
end

