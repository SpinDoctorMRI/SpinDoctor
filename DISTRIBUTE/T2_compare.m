clear all;
exno=1;
markers={'x','o','s','v','t'};
myfolder='./';
for i=1:3
    if i==1
        load([myfolder,'T2_20e3_20e3_20e3_mp_dt1000']);
    end;
    if i==1
%         load([myfolder,'T2_40e3_40e3_40e3']);
        load([myfolder,'T2_20e3_20e3_20e3_new1']);
    end;
    if i==2
        load([myfolder,'T2_20e3_20e3_20e3_new']);
    end;
%     if i==4
%         load([myfolder,'T2_20e3_20e3_20e3']);
%     end;
%     if i==5
%         load([myfolder,'T2_10e3_10e3_10e3']);
%     end;
%     if i==6
%         load([myfolder,'T2_40e3_20e3_80e3']);
%     end;
    disp([num2str(i),') ',num2str(T2_cmpts)]);
%     Sig_free_new = reshape(Sig_free, size(experi_btpde.bvalues));
%     PLOT_SIGNAL(experi_btpde.bvalues(exno,:),MF_allcmpts(exno,:),Sig_free_new(exno,:),ADC_allcmpts_S0(exno),ADC_allcmpts(exno))
    NormalizedMF=MF_allcmpts./(ones(size(MF_allcmpts))*MF_allcmpts(exno,1));
    semilogy(experi_btpde.bvalues(exno,:), real(NormalizedMF(exno,:)),markers{i},'MarkerSize',10);
%     kk=1;
%     if (i==1)
%         TE = experi_btpde.sdeltavec(exno) + experi_btpde.bdeltavec(exno);
%         kk = exp(TE/T2_cmpts(1));
%         [TE, kk]
%     end;
%     semilogy(experi_btpde.bvalues(exno,:), kk*real(MF_allcmpts(exno,:)));
    hold all;
    set(gca,'FontSize',15)
    xlabel('b (s/mm^2)');
    ylabel('S(b)/S(0)');
    
% %     ylim([2.5, 3.6])
%     format short e
%     ADC_allcmpts
end;


