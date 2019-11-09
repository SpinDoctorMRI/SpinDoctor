function [EIG_length_cmpts]=MF_EIG_TO_LENGTH(EIG_value_cmpts,DIFF_cmpts)

% convert the computed eigenvalues into a length scale
% 
% Input:
%     1. EIG_value_cmpts
%     2. DIFF_cmpts
% 
% Output:
%     EIG_length_cmpts

cvec = [1.8412,3.0542];
rmax = 0;
hmax = 0;

Ncmpt = length(DIFF_cmpts);

for icmpt = 1:Ncmpt
    neig = length(EIG_value_cmpts{icmpt});
    EIG_length_cmpts{icmpt} = inf*ones(neig,1);
    ii = find(EIG_value_cmpts{icmpt}>1e-16);    
    lvec = EIG_value_cmpts{icmpt}(ii)/DIFF_cmpts(icmpt);
    rvec = sqrt(cvec(1)^2./lvec);
    hvec = sqrt(1^2*pi^2./lvec);
    EIG_length_cmpts{icmpt}(ii) = hvec;
end

