% function Tm = get_Tm_camino(seq)
%GET_TM_CAMINO computes the temporal matrix for a camino sequence for generalised Mitra formula.
%   Moutal, N., Maximov, I.I. and Grebenkov, D.S., 2019. 
%   Probing surface-to-volume ratio of an anisotropic medium by diffusion NMR with general gradient encoding. 
%   IEEE Transactions on Medical Imaging, 38(11), pp.2507-2522.

% Note some problem need to double check it is correct.

seq = read_scheme('camino_sequences/test_scheme.scheme');
seq  = seq{1};


%%
% clear seq;
% K = 101;
% dt = 1e-2;
% t = dt*(0:(K-1));
% g   = [t;t;t];
% seq = SequenceCamino(K,dt,g,'test');

%%
g = seq.g;
K = seq.K; dt = seq.dt;
T = seq.echotime;
g = (g(:,1:(K-1)) + g(:,2:end))/2;


%%

y = tensorprod(g,g');

t = dt*(0:K-1);

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

Int = zeros(3,3);

for i = 1:K-1
    for j=1:K-1
        Int = Int + squeeze(W(i,j).*y(:,i,j,:));
    end
end


%%
gamma = 0.2675222005 *1e-6;
B = seq.get_B_tensor(gamma);
b = trace(B);
Tm = -(gamma^2*T)/(2*b)*Int;