function [D,root_mean_sq_err] = fit_angular_diffusion_field(S,b)
%FIT_ANGULAR_DIFFUSION_FIELD finds a function D : S^2 -> R which approximatley satisfies
% -D(u)b = ln(S(u,b))
% Input:
%       S: volume weighted signal
%       b: b values
% Outpur:
%       D: diffusivity field


% d = sum(b.^2)*sum(b.^4) - sum(b.^3)^2;
% 
% A = 1/d * [sum(b.^4), - sum(b.^3) ;-sum(b.^3) , sum(b.^2)];
% V = [(- b'*log(S))',( - (b.^2)'*log(S))'];
% D = A*V';
% D =D';
% % 
D = - b'*log(S)./sum(b.^2);
root_mean_sq_err = vecnorm(b*D +log(S),2,1)/sqrt(length(b));
