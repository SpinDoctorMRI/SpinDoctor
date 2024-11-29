% function Tm = get_Tm_PGSE(seq,g,ug,b)
%GET_TM_CAMINO computes the temporal matrix for a PGSE sequence for generalised Mitra formula.
%   Moutal, N., Maximov, I.I. and Grebenkov, D.S., 2019. 
%   Probing surface-to-volume ratio of an anisotropic medium by diffusion NMR with general gradient encoding. 
%   IEEE Transactions on Medical Imaging, 38(11), pp.2507-2522.

seq = PGSE(100,2000); g = 100; ug = [1;0;0];

d = seq.delta; D = seq.Delta;
T = seq.TE;

t= [0 d d+D T];

% t = linspace(0,T,100);
% y = seq.call(t);
% plot(t,y)


K = 4;


W =  zeros(K-1,K-1);

m = 3;

for i = 2:K
    for j=2:i
        W(i-1,j-1) = 4*T^(-m)/((m+2)*(m+4)) * (-abs(t(i) - t(j))^(m/2 + 2) + ...
        abs(t(i-1) - t(j))^(m/2 + 2) + ...
        abs(t(i) - t(j-1))^(m/2 + 2) + ...
        -abs(t(i-1) - t(j-1))^(m/2 + 2)  );
        W(j-1,i-1) = W(i-1,j-1);
    end
end

Int_1d = (W(1,1) - W(1,3) - W(3,1) + W(3,3))*g^2;


gamma = 0.2675222005 *1e-6;
Tm = -(gamma^2*T)/(2*b)*Int_1d *(ug*ug');
