function [Deff_eig_tensor_cmpts,Deff_eig_tensor_allcmpts] ...
    = MF_DTENSOR(VOL_cmpts,IC_cmpts,EIG_value_cmpts,EIG_proj_cmpts,DIFF_cmpts,MF_JN_cmpts)

% compute the effective diffusion tensor
%
% Input:
%     1. VOL_cmpts
%     2. IC_cmpts
%     3. EIG_value_cmpts
%     4. EIG_proj_cmpts
%     5. DIFF_cmpts
%     6. MF_JN_cmpts
%
% Output:
%     1. Deff_eig_tensor_cmpts
%     2. Deff_eig_tensor_allcmpts


disp(['DEFF MF: eigenfunction decomposition']);

Ncmpt = length(VOL_cmpts);

VOL_allcmpts = sum(VOL_cmpts);
VOL_frac = VOL_cmpts/VOL_allcmpts;

nexperi = size(MF_JN_cmpts{1},1);
elapsed_time = zeros(Ncmpt, nexperi);


Deff_eig_tensor_cmpts = zeros(Ncmpt,nexperi,3,3);

for icmpt = 1:Ncmpt  
    neig = length(EIG_value_cmpts{icmpt});
   
    for iexperi = 1:nexperi
        Deff_eig_tensor{icmpt} = zeros(neig,3,3);
        disp(['      Experiment ',num2str(iexperi)]);
        for ieig = 1:neig
            [dtmp]=[MF_JN_cmpts{icmpt}(iexperi,ieig)*DIFF_cmpts(icmpt)];      
            vvec = squeeze(EIG_proj_cmpts{icmpt}(1:3,1,ieig));        
            Deff_eig_tensor{icmpt}(ieig,:,:) = dtmp*vvec*vvec';
        end
        Deff_eig_tensor_cmpts(icmpt,iexperi,:,:) = squeeze(sum(Deff_eig_tensor{icmpt},1));        
    end
end

Deff_eig_tensor_allcmpts = zeros(nexperi,3,3);
for iexperi = 1:nexperi
    for icmpt = 1:Ncmpt
        Deff_eig_tensor_allcmpts(iexperi,:,:) = squeeze(Deff_eig_tensor_allcmpts(iexperi,:,:)) + ...
            (IC_cmpts(icmpt)*VOL_frac(icmpt)*squeeze(Deff_eig_tensor_cmpts(icmpt,iexperi,:,:)))/sum((IC_cmpts.*VOL_frac)');
    end
end