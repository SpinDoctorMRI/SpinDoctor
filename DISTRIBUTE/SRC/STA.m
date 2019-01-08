function [ADC_STA,ADC_STA_allcmpts] = STA(experiment,DIFF_cmpts,...
            VOL,SAu,IC_cmpts)

			
sdeltavec = experiment.sdeltavec;
bdeltavec = experiment.bdeltavec;
seqvec = experiment.seqvec;
npervec = experiment.npervec;
nexperi = length(sdeltavec);
Ncmpt = length(DIFF_cmpts);

VOL_allcmpts = 0;

for icmpt = 1:Ncmpt
    VOL_allcmpts  = VOL_allcmpts + VOL(icmpt);
end

for icmpt = 1:Ncmpt
    VOL_frac(icmpt) = VOL(icmpt)/VOL_allcmpts;
end

ADC_STA = zeros(Ncmpt,nexperi);

for iexperi = 1:nexperi
    for icmpt = 1:Ncmpt
        [ADC_STA(icmpt,iexperi)] = deff_sta(DIFF_cmpts(icmpt),...
            VOL(icmpt),SAu(icmpt),sdeltavec(iexperi),bdeltavec(iexperi),...
            seqvec(iexperi),npervec(iexperi));     
    end
end	
ADC_STA_allcmpts = nan*ones(nexperi,1);
for iexperi = 1:nexperi
    ADC_STA_allcmpts(iexperi,1) = sum((IC_cmpts.*VOL_frac)'.*ADC_STA(:,iexperi))./sum((IC_cmpts.*VOL_frac)');
end