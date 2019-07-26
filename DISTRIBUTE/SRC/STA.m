function [ADC_STA_cmpts,ADC_STA_allcmpts] = STA(experiment,DIFF_cmpts,...
    VOL_cmpts,SAu_cmpts,IC_cmpts)

% compute the ADC in the short diffusion time regime
% 
% Input:
%     1. experiment is a structure with 8 elements:
%         ngdir_total 
%         gdir        
%         sdeltavec   
%         bdeltavec   
%         seqvec      
%         npervec    
%         rtol       
%         atol        
%     2. DIFF_cmpts
%     3. VOL_cmpts
%     4. SAu_cmpts
%     5. IC_cmpts
% 
% Output: 
%     1. ADC_STA_cmpts
%     2. ADC_STA_allcmpts

sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;
nexperi = length(sdeltavec);
Ncmpt = length(DIFF_cmpts);

VOL_allcmpts = 0;

for icmpt = 1:Ncmpt
    VOL_allcmpts  = VOL_allcmpts + VOL_cmpts(icmpt);
end

for icmpt = 1:Ncmpt
    VF_cmpts(icmpt) = VOL_cmpts(icmpt)/VOL_allcmpts;
end

ADC_STA_cmpts = zeros(Ncmpt,nexperi);

for iexperi = 1:nexperi
    for icmpt = 1:Ncmpt
        [ADC_STA_cmpts(icmpt,iexperi)] = deff_sta(DIFF_cmpts(icmpt),...
            VOL_cmpts(icmpt),SAu_cmpts(icmpt),sdeltavec(iexperi),bdeltavec(iexperi),...
            seqvec(iexperi),npervec(iexperi));
    end
end
ADC_STA_allcmpts = nan*ones(nexperi,1);
for iexperi = 1:nexperi
    ADC_STA_allcmpts(iexperi,1) = sum((IC_cmpts.*VF_cmpts)'.*ADC_STA_cmpts(:,iexperi))./sum((IC_cmpts.*VF_cmpts)');
end