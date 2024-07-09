 L_max = 8;
% figure;
% t = tiledlayout(2,4);
for R0 = 8
r = unitsphere(64);
l = 0:2:L_max;
I = zeros(length(l),length(setup.gradient.directions));
p = zeros(L_max + 1,2*L_max + 1);
for j = 1:length(l)
    I(j,:) = get_I_l(D,R0,l(j),setup);
    a_l = harmonic_analysis_laplace(I(j,:),u,L_max);
    p(l+1,:) = (-1i)^(l(j)/2).*a_l(l+1,:);
end

Y = zeros(L_max + 1,2*L_max + 1,length(r));
for j=1:length(l)
    for m = -l(j):l(j)
        Y(1+l(j),L_max + 1 + m,:) = harmonic(r',l(j),m) ;
    end
end

V = norm(reshape(p(3:2:end,:),[],1)).^2/(9*abs(p(1,L_max+1)).^2);
Integral = real(p(1,L_max+1));
Y = reshape(Y,[],length(r));
p = reshape(p,1,[]);

P = real(p*Y)./(sqrt(4*pi)*Integral);
% nexttile;
% plot_scalar_field(r,P);
% title(sprintf('$R_0 = %d$',R0),'Interpreter','latex');
end
% title(t,sprintf('%s',cellname),'Interpreter','none');
