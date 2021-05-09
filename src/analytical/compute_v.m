function [v, dv] = compute_v(alpha, n, x, params, bc)
%COMPUTE_V Compute the eigenfunction v and its derivative.

r = params.r;
D = params.D;
d = params.d;

m = length(D);

v = ones(size(x));
dv = zeros(size(x));

if alpha > 0
    if nargin < 5
        bc = compute_bc(alpha, n, params);
    end

    r1(2:m+1) = r(1:m);
    r1(1) = -1e-8;

    for i = 1:m
        inds = r1(i) < x & x <= r1(i + 1);
        [J, Y, dJ, dY] = compute_JY(alpha * x(inds) / sqrt(D(i)), n, d);
        if abs(bc(i, 2)) > 1e-12
            v(inds) = bc(i, 1) * J + bc(i, 2) * Y;
            dv(inds) = alpha / sqrt(D(i)) * (bc(i, 1) * dJ + bc(i, 2) * dY);
        else
            v(inds) = bc(i, 1) * J;
            dv(inds) = alpha / sqrt(D(i)) * bc(i, 1) * dJ;
        end
    end
end
