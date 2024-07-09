function [f] = expmv_old(A, b, s, m, mu, tol)
%EXPMV   Matrix exponential times vector or matrix.
% MODIFIED BASED ON Higham's ORIGNAL VERSION.
%   [F,S,M,MV,MVD] = EXPMV(t,A,B,[],PREC) computes EXPM(t*A)*B without
%   explicitly forming EXPM(t*A). PREC is the required accuracy, 'double',
%   'single' or 'half', and defaults to underlyingType(A).
%   A total of MV products with A or A^* are used, of which MVD are
%   for norm estimation.
%   The full syntax is
%     [f,s,m,mv,mvd] = expmv(t,A,b,M,prec,shift,bal,full_term,prnt).
%   unA = 1 if the alpha_p were used instead of norm(A).
%   If repeated invocation of EXPMV is required for several values of t
%   or B, it is recommended to provide M as an external parameter as
%   M = SELECT_TAYLOR_DEGREE(A,b,m_max,p_max,prec,shift,bal,true).
%   This also allows choosing different m_max and p_max.

%   Reference: A. H. Al-Mohy and N. J. Higham. Computing the action of the
%   matrix exponential, with an application to exponential integrators.
%   SIAM J. Sci. Comput., 33(2):488--511, 2011.  Algorithm 3.2.

%   Awad H. Al-Mohy and Nicholas J. Higham, November 9, 2010.

full_term = false;

eta = exp(mu/s);
f = b;
for ii = 1:s
    c1 = norm(b, inf);
    for k = 1:m
        b = (1/(s*k))*(A*b);
        f =  f + b;
        c2 = norm(b,inf);
        if ~full_term
            if c1 + c2 <= tol*norm(f,inf)
                break
            end
            c1 = c2;
        end
    end
    f = eta*f; b = f;
end
