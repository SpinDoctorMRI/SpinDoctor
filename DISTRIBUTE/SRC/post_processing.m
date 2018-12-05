function [ADC_allcmpts,ADC_allcmpts_polydeg,ADC_allcmpts_S0,...
    ADC,ADC_polydeg,ADC_S0,MF_allcmpts,M0_allcmpts,S0_allcmpts,...
    MF,M0,S0]...
    = post_processing(MT,bvalues)

nexperi = length(MT);
nb = length(MT{1});
Ncmpt = length(MT{1}{1});
for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);  
    nb = length(bvec);
    for ib = 1:nb
        for icmpt = 1:Ncmpt
            MF(icmpt,iexperi,ib) = MT{iexperi}{ib}{icmpt}(end);
            M0(icmpt,iexperi,ib) = MT{iexperi}{ib}{icmpt}(1);
        end
        MF_allcmpts(iexperi,ib) = 0;
        for icmpt = 1:Ncmpt
            MF_allcmpts(iexperi,ib) = MF_allcmpts(iexperi,ib) + MF(icmpt,iexperi,ib);
        end
        M0_allcmpts(iexperi,ib) = 0;
        for icmpt = 1:Ncmpt
            M0_allcmpts(iexperi,ib) = M0_allcmpts(iexperi,ib) + M0(icmpt,iexperi,ib);
        end
    end   
    ib0 = find(abs(bvec)<=1e-16);
    ibn0 = find(abs(bvec)>1e-16);
    if (length(ib0) >= 1)
        for icmpt = 1:Ncmpt
            S0(icmpt,iexperi) = MF(icmpt,iexperi,ib0(1));
        end
        S0_allcmpts(iexperi) = MF_allcmpts(iexperi,ib0(1));
    else
        S0(1:Ncmpt,iexperi) = nan;
        S0_allcmpts(iexperi) = nan;
    end
end
ADC = nan*ones(Ncmpt,nexperi);
ADC_allcmpts = nan*ones(nexperi,1);
ADC_polydeg = nan*ones(Ncmpt,nexperi);
ADC_allcmpts_polydeg = nan*ones(nexperi,1);
for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);
    if (length(bvec) >= 2)       
        bmin = bvec(1);
        bmax = bvec(end);        
        for icmpt = 1:Ncmpt
            data1d = real(squeeze(MF(icmpt,iexperi,:)))';
            
            [fit_poly,ADC01d,KUR1d,KUR01d,S01d,Cfit1d,errfit,ndeg,ADC0_err1d,KUR_err1d] ...
                = process_signal_POLY(data1d,bvec,bmin,bmax);
            ADC(icmpt,iexperi) = ADC01d;
            ADC_polydeg(icmpt,iexperi) = ndeg;
            ADC_S0(icmpt,iexperi) = S01d;
        end
        data1d = real(MF_allcmpts(iexperi,:));
        [fit_poly,ADC01d,KUR1d,KUR01d,S01d,Cfit1d,errfit,ndeg,ADC0_err1d,KUR_err1d] ...
            = process_signal_POLY(data1d,bvec,bmin,bmax);
        ADC_allcmpts(iexperi,1) = ADC01d;
        ADC_allcmpts_polydeg(iexperi,1) = ndeg;
        ADC_allcmpts_S0(iexperi,1) = S01d;
    end
end






