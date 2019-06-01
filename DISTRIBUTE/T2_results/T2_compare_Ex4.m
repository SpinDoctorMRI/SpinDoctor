clear all; path(pathdef); clc;
exno=1;
markers={'--o','--p','--x','--v','--s'};
myfolder='./';

for kk=1:2
    seq='PGSE';
    if kk==2
            seq='dPGSE';
    end;
    for i=1:4
        if i==1
            load([myfolder,'T2_Inf_Inf_Inf_',seq]);
        end;    
        if i==2
            load([myfolder,'T2_40e3_40e3_40e3_',seq]);
        end;
        if i==3
            load([myfolder,'T2_20e3_20e3_20e3_',seq]);
        end;
        if i==4
            load([myfolder,'T2_40e3_20e3_80e3_',seq]);
        end;
        disp([num2str(i),') ',num2str(T2_cmpts)]);
        NormalizedMF=MF_allcmpts./(ones(size(MF_allcmpts))*MF_allcmpts(exno,1));
        subplot(2,2, 2*(kk-1) + 2);
        semilogy(experi_btpde.bvalues(exno,:), real(NormalizedMF(exno,:)),markers{i},'MarkerSize',10);
        box off; hold all;
        set(gca,'FontSize',15)        
        ylabel([seq,', S(b)/S(0)']);
        ylim([5e-2, 1]);
        yticks([5e-2, 1e-1, 5e-1, 1e0])
  
        subplot(2,2, 2*(kk-1)+1);
        semilogy(experi_btpde.bvalues(exno,:), real(MF_allcmpts(exno,:)), markers{i},'MarkerSize',10);
        ylabel([seq,', S(b)']);
        set(gca,'FontSize',15);
        ylim([1e1, 1e4]);
        yticks([1e1, 1e2, 1e3, 1e4])
        box off; hold all;
    end;
end;

for k=1:2
    subplot(2,2,k);
    xlabel('b (s/mm^2)');
end;

legend({'T2=[\infty,\infty,\infty]','T2=[40,40,40] ms', 'T2=[20,20,20] ms', 'T2=[40,20,80] ms'}, 'Location', 'southoutside', 'NumColumns',4,'FontSize',14);
legend boxoff