function [cost, A_shift, s, m, mu, tol] = expmv_cost(A, b, M)
%EXPMV   EXPMV computational cost.
% MODIFIED BASED ON Higham's ORIGNAL VERSION.

%   Reference: A. H. Al-Mohy and N. J. Higham. Computing the action of the
%   matrix exponential, with an application to exponential integrators.
%   SIAM J. Sci. Comput., 33(2):488--511, 2011.  Algorithm 3.2.

%   Awad H. Al-Mohy and Nicholas J. Higham, November 9, 2010.

prec = underlyingType(A);
n = length(A);
t = 1;

switch prec
    case 'double', tol = 2^(-53);
    case 'single', tol = 2^(-24);
    case 'half',   tol = 2^(-10);
end

% shift
mu = gather(trace(A)/n); 
ss = size(A);
% Indices of the main diagonal
index = 1:ss(1)+1:ss(1)*ss(2);
A_shift = A;
A_shift(index) = A_shift(index) - mu;

if nargin < 3 || isempty(M)
    tt = 1;
    [M,mvd,alpha] = select_taylor_degree(A_shift, b, [], [], prec, false, false);
    mv = mvd;
else
    tt = t; mv = 0; mvd = 0;
end

s = 1;
if t == 0
    m = 0;
else
    [m_max, p] = size(M);
    U = diag(1:m_max);
    C = (ceil(abs(tt)*M))'*U;
    C(C == 0) = inf;
    if p > 1
        % cost is the overall cost
        [cost, m] = min(min(C));
    else
        % when C is one column. Happens if p_max = 2
        [cost, m] = min(C);
    end
    if cost == inf; cost = 0; end
    s = max(cost/m,1);
end
