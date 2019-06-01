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
            load([myfolder,'T2_40e3_20e3_20e3_',seq]);
        end;
        disp([num2str(i),') ',num2str(T2_cmpts)]);
        NormalizedMF=MF_allcmpts./(ones(size(MF_allcmpts))*MF_allcmpts(exno,1));
        subplot(2,2, 2*(kk-1) + 2);
        semilogy(experi_btpde.bvalues(exno,:), real(NormalizedMF(exno,:)),markers{i},'MarkerSize',10);
        box off; hold all;
        ylabel([seq,', S(b)/S(0)']);

        subplot(2,2, 2*(kk-1)+1);
        semilogy(experi_btpde.bvalues(exno,:), real(MF_allcmpts(exno,:)), markers{i},'MarkerSize',10);
        ylabel([seq,', S(b)']);

        box off; hold all;
    end;
end;

for k=1:4
    subplot(2,2,k);
    set(gca,'FontSize',15)
    xlabel('b (s/mm^2)');
end;

legend({'T2=[\infty, \infty]','T2=[40,40] ms', 'T2=[20, 20] ms', 'T2=[40, 20] ms'}, 'Location', 'southoutside', 'NumColumns',4);
