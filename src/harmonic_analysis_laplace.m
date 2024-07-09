function [a,C,ODF] = harmonic_analysis_laplace(S,u,L_max)
%HARMONIC_ANALYSIS_LAPLACE computes the SHT of a signal sampled at points u
%on the unit sphere using a Laplace-Betrami regularisation operator.
% 
% Descoteaux, M., Angelino, E., Fitzgibbons, S. and Deriche, R. (2007), 
% Regularized, fast, and robust analytical Q-ball imaging. 
% Magn. Reson. Med., 58: 497-510. https://doi.org/10.1002/mrm.21277
S = real(S); 
R = (L_max+1)*(L_max+2)/2;
 N = length(S);
 Yjs = cell(R,1);
 B = zeros(N,R);
 L = zeros(R,1);
 i = zeros(R,1);
 k = zeros(R,1);
 for l = 0:2:L_max
     for m = -l:l
         j = J(l,m);
         i(j) = l;
         k(j) = m;   
         if m <=0
            Yjs{j} = @(p) sqrt(2)*real(harmonic(p,l,m));
         else
            Yjs{j} = @(p) sqrt(2)*imag(harmonic(p,l,m));
         end
         L(j) = l^2*(l+1)^2;
     end
 end
L = diag(L);
S = reshape(S,[],1);

for r = 1:R
    B(:,r) = Yjs{r}(u');
end
lambda = 0.006;
M = (B'*B+lambda*L);
C = M \ (B'*S);

a = zeros(L_max+1,2*L_max+1);
for r = 1:R
 a(1+i(r),L_max +1+ k(r))= C(r);
end

Z2 = 0;ODF=0;
for r = 1:R
    l=i(r); j = k(r);
    Pl= (-1)^(l/2)*prod(1:2:(l-1))/prod(2:2:l);
    ODF = ODF+2*pi*C(r)*Pl*Yjs{r}(u');
    Z2 = Z2 + 4*pi^2*abs(C(r))^2;
end

ODF = ODF/sqrt(Z2);


end



function j = J(l,m)
    j = (l^2 + l+ 2)/2 +m;
end