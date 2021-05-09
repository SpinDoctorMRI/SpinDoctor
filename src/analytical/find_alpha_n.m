function alpha = find_alpha_n(alpha_max, n, params, dalpha)
%FIND_ALPHA_N Find all radial eigenvalues for a given angular mode n.
%   This function sweeps over the interval [0, alpha_max].
%
%   alpha_max: [1 x 1]
%   n
%   params
%   dalpha


alpha = [];

% Lower bound (do not include zero)
if n < 20
    a = 0.01;
elseif n < 40
    a = 1;
elseif n < 60
    a = 10;
else
    a = 100;
end

F = alpha_func(a, n, params);
assert(~isnan(F), "Too small initial value a");

% Find zeros of Fn(alpha)
k = 1;
while a < alpha_max
    F_old = F;
    a_old = a;
    
    % Keep advancing alpha until F(alpha) changes sign or max alpha is attained
    while F_old * F >= 0 && a < alpha_max
        a = a + dalpha;
        F = alpha_func(a, n, params);
    end
    if a < alpha_max
        % Find alpha by dichotomy between previous alpha new upper bound
        alpha(k) = find_alpha_n_interval(a_old, a, 0, n, params);
    end
    
    k = k + 1;
end
