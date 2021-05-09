function f = alpha_func(alpha, n, params)
%ALPHA_FUNC Compute the quantity Fn(alpha), for which the roots are the radial eigenvalues.
%
%   alpha
%   n
%   params
%
%   f

r = params.r;
D = params.D;
W = params.W;
d = params.d;

m = length(D);
bc1 = zeros(size(alpha));
bc2 = zeros(size(alpha));

for i = 1:length(alpha)
    bc = compute_bc(alpha(i), n, params);
    bc1(i) = bc(m, 1);
    bc2(i) = bc(m, 2);
end

rD = r(m) / sqrt(D(m));

if alpha == 0
    error('alpha zero');
end

switch d
    case 2
        f1 = sqrt(D(m)) * alpha .* (besselj(n - 1, alpha * rD) ...
            - besselj(n + 1, alpha * rD)) / 2 ...
            + W(m) * besselj(n, alpha * rD);
        f2 = sqrt(D(m)) * alpha .* (bessely(n - 1, alpha * rD) ...
            - bessely(n + 1, alpha * rD)) / 2 ...
            + W(m) * bessely(n, alpha * rD);
    case 3
        f1 = sqrt(D(m)) * alpha .* (besseljs(n - 1, alpha * rD) ...
            - besseljs(n + 1, alpha * rD) ...
            - besseljs(n, alpha * rD) ./ (alpha * rD) ) / 2 ...
            + W(m) * besseljs(n, alpha * rD);
        f2 = sqrt(D(m)) * alpha .* (besselys(n - 1, alpha * rD) ...
            - besselys(n + 1, alpha * rD) ...
            - besselys(n, alpha * rD) ./ (alpha * rD) ) / 2 ...
            + W(m) * besselys(n, alpha * rD);
end

f = bc1 .* f1 + bc2 .* f2;
