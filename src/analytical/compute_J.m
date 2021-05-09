function J = compute_J(alpha, params)
%COMPUTE_J Compute the integral quantity J.

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
