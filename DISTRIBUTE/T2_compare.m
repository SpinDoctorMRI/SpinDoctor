clear all;
exno=1;

myfolder='./';
for i=1:2
    if i==1
        load([myfolder,'T2_Inf_Inf_Inf_new']);
    end;
    if i==2
        load([myfolder,'T2_20e3_20e3_20e3_new']);
    end;
    if i==3
        load([myfolder,'T2_40e3_40e3_40e3']);
    end;
    if i==4
        load([myfolder,'T2_20e3_20e3_20e3']);
    end;
    if i==5
        load([myfolder,'T2_10e3_10e3_10e3']);
    end;
    if i==6
        load([myfolder,'T2_40e3_20e3_80e3']);
    end;
    disp([num2str(i),') ',num2str(T2_cmpts)]);
    Sig_free_new = reshape(Sig_free, size(experi_btpde.bvalues));
%     PLOT_SIGNAL(experi_btpde.bvalues(exno,:),MF_allcmpts(exno,:),Sig_free_new(exno,:),ADC_allcmpts_S0(exno),ADC_allcmpts(exno))
%     NormalizedMF=MF_allcmpts./(ones(size(MF_allcmpts))*MF_allcmpts(exno,1));
%     semilogy(experi_btpde.bvalues(exno,:), real(NormalizedMF(exno,:)));
%     kk=1;
%     if (i==1)
%         TE = experi_btpde.sdeltavec(exno) + experi_btpde.bdeltavec(exno);
%         kk = exp(-TE/20000);
%         [TE, kk]
%     end;
    semilogy(experi_btpde.bvalues(exno,:), real(MF_allcmpts(exno,:)));
    hold all;
%     ylim([2.5, 3.6])
    format short e
    ADC_allcmpts
end;


