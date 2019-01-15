function M=mass_matrixP1_2D(elements,areas,coeffs)
% Copyright (c) 2013, Talal Rahman, Jan Valdman

%coeffs can be only P0 (elementwise constant) function 
%represented by a collumn vector with size(elements,1) entries
%if coeffs is not provided then coeffs=1 is assumed globally
%Note: P1 coeffs needs a higher integration rule (not implemented yet)

%KRON   Kronecker tensor product.
%   KRON(X,Y) is the Kronecker tensor product of X and Y.
%   The result is a large matrix formed by taking all possible
%   products between the elements of X and those of Y. For
%   example, if X is 2 by 3, then KRON(X,Y) is
%
%      [ X(1,1)*Y  X(1,2)*Y  X(1,3)*Y
%        X(2,1)*Y  X(2,2)*Y  X(2,3)*Y ]
%
%   If either X or Y is sparse, only nonzero elements are multiplied
%   in the computation, and the result is sparse.
%


%   SPARSE Create sparse matrix.
%   S = SPARSE(X) converts a sparse or full matrix to sparse form by squeezing
%   out any zero elements.
%
%   S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to generate an
%   m-by-n sparse matrix such that S(i(k),j(k)) = s(k), with space
%   allocated for nzmax nonzeros.  Vectors i, j, and s are all the same
%   length.  Any elements of s that are zero are ignored, along with the
%   corresponding values of i and j.  Any elements of s that have duplicate
%   values of i and j are added together.  The argument s and one of the
%   arguments i or j may be scalars, in which case the scalars are expanded
%   so that the first three arguments all have the same length.
%
%   S = SPARSE(i,j,s,m,n) where nzmax = length(s).
%
%   S = SPARSE(i,j,s) where m = max(i) , n = max(j) and  nzmax = length(s).
%
X=kron(ones(1,3),elements); Y=kron(elements,ones(1,3)); 

if (nargin<3)
    Zmassmatrix=kron(areas,reshape((ones(3)+eye(3))/12,1,9)); 
else
    if numel(coeffs)==size(elements,1) %P0 coefficients
        Zmassmatrix=kron(areas.*coeffs,reshape((ones(3)+eye(3))/12,1,9)); 
    else %P1 coefficients
        M1=[6 2 2; 2 2 1; 2 1 2]/60;
        M2=M1([3,1,2],[3,1,2]);
        M3=M2([3,1,2],[3,1,2]);
            
        Zmassmatrix=kron(areas.*coeffs(elements(:,1)),reshape(M1,1,9)) ...
                   +kron(areas.*coeffs(elements(:,2)),reshape(M2,1,9)) ...
                   +kron(areas.*coeffs(elements(:,3)),reshape(M3,1,9));
    end
        
end

M=sparse(X,Y,Zmassmatrix); %M(X(k),Y(k)) = Z(k)
