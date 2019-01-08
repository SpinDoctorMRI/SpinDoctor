function [Sig_free,ADC_free_allcmpts] = ADCFREE(bvalues,DIFF_cmpts,VOL,IC_cmpts)

nexperi = length(bvalues);
Ncmpt = length(DIFF_cmpts);

VOL_allcmpts = 0;

for icmpt = 1:Ncmpt
    VOL_allcmpts  = VOL_allcmpts + VOL(icmpt);
end

for icmpt = 1:Ncmpt
    VOL_frac(icmpt) = VOL(icmpt)/VOL_allcmpts;
end

Sig_free = zeros(size(bvalues(:)));

for icmpt = 1:Ncmpt
    Sig_free = Sig_free+IC_cmpts(1,icmpt)*VOL(icmpt)*exp(-DIFF_cmpts(icmpt)*bvalues(:));
end

ADC_free_allcmpts = sum((IC_cmpts.*VOL_frac)'.*DIFF_cmpts')./sum((IC_cmpts.*VOL_frac)');

