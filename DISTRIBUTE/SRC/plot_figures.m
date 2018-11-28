function plot_figures(UG,bvalues,difftime,...
    DIFF_cmpts,IC_cmpts,Ncmpt,ADC_allcmpts,ADC_allcmpts_S0,...
    Deff_STA_allcmpts,ADC,ADC_S0,...
    Deff_STA,MF_allcmpts,MF,VOL,SAu,VOL_frac)

mydefinitions;

ncolor = length(colorvec_cell);
 
Sig_free = zeros(size(bvalues));

nexperi = length(difftime);

for iexperi = 1:nexperi
    vol = 0;
    for icmpt = 1:Ncmpt
        Sig_free(iexperi,:) = Sig_free(iexperi,:)+IC_cmpts(1,icmpt)*VOL_frac(icmpt)*exp(-DIFF_cmpts(icmpt)*bvalues(iexperi,:));
        vol = vol+IC_cmpts(1,icmpt)*VOL_frac(icmpt);
    end
    Sig_free(iexperi,:) = Sig_free(iexperi,:)/vol;
end

figure; hold on
iplot = 0;
for icmpt = 1:Ncmpt
    for iexperi = 1:nexperi
        bvec = bvalues(iexperi,:);
        yvec = real(squeeze(MF(icmpt,iexperi,:)));
        h = plot(bvec, log10(yvec), ...
            [colorvec_cell{mod(icmpt-1,ncolor)+1},markervec_cell{iexperi}]);
        set(h,'MarkerSize', 10, 'LineWidth',1);
        iplot = iplot + 1;
        legend_vec{iplot} = ['Cmpt ',num2str(icmpt),', Exp ',mynum2str(iexperi)];
    end    
end
for iexperi = 1:nexperi
    for icmpt = 1:Ncmpt
        yvec = ADC_S0(icmpt,iexperi)*exp(-ADC(icmpt,iexperi)*(bvec));
        h = plot(bvec, log10(yvec), ...
            [colorvec_cell{mod(icmpt-1,ncolor)+1},'-']);
        set(h,'MarkerSize', 10, 'LineWidth',1);
    end
end
legend(legend_vec{1:iplot},'Location','NorthEastOutside');
set(legend, 'FontSize',10)
set(gca, 'FontSize',10)
xlabel('bvalue')
ylabel('Sig')
title(['UG = [',mynum2str(UG),']'])



figure; hold on
iplot = 0;

for iexperi = 1:nexperi   
    yvec = real(MF_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{1},markervec_cell{iexperi}]);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = [' Exp ',mynum2str(iexperi)];

    yvec = ADC_allcmpts_S0(iexperi)*Sig_free(iexperi,:);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{2},markervec_cell{iexperi},'-']);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = ['free diffusion'];
end 

for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);
    yvec = ADC_allcmpts_S0(iexperi)*exp(-ADC_allcmpts(iexperi)*bvec);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{1},'-']);
    set(h,'MarkerSize', 10, 'LineWidth',1);
end
legend(legend_vec{1:iplot},'Location','NorthEastOutside');
set(legend, 'FontSize',10)
set(gca, 'FontSize',10)
xlabel('bvalue')
ylabel('Sig')
title(['UG = [',mynum2str(UG),']'])

figure; hold on
iplot = 0;
xvec = SAu./VOL;
for iexperi = 1:nexperi
    yvec = squeeze(ADC(:,iexperi));
    h = plot(xvec,yvec,[colorvec_cell{mod(iexperi-1,ncolor)+1},markervec_cell{1}]);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = ['Exp ',num2str(iexperi), ', ', 'SIMUL'];
    
    yvec = squeeze(Deff_STA(:,iexperi));
    h = plot(xvec,yvec,[colorvec_cell{mod(iexperi-1,ncolor)+1},markervec_cell{2}]);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = ['Exp ',num2str(iexperi), ', ', 'Deff STA'];
end
legend(legend_vec{1:iplot},'Location','NorthEastOutside');
set(legend, 'FontSize',10)
set(gca, 'FontSize',10)
xlabel('SoV')
ylabel('D^{eff}')
title(['UG = [',mynum2str(UG),']'])
set(gca,'xlim',[0,max(xvec)*1.2]);
set(gca,'ylim',[0,max(DIFF_cmpts)*1.1]);



figure; hold on

iplot = 0;
xvec = sqrt(difftime/1000);
yvec = ADC_allcmpts(:);
h = plot(xvec,yvec,[colorvec_cell{1},'d']);
set(h,'MarkerSize', 10, 'LineWidth',1);
iplot = iplot + 1;
legend_vec{iplot} = ['ADC Allcmpts'];

yvec = Deff_STA_allcmpts;
h = plot(xvec,yvec,[colorvec_cell{2},'d']);
set(h,'MarkerSize', 10, 'LineWidth',1);
iplot = iplot + 1;
legend_vec{iplot} = ['Deff STA Allcmpts'];

legend(legend_vec{1:iplot},'Location','NorthEastOutside');
set(legend, 'FontSize',10)
set(gca, 'FontSize',10)
xlabel('sqrt (Diffusion time (ms))')
ylabel('D^{eff}')
title(['UG = [',mynum2str(UG),']'])
set(gca,'xlim',[0,max(xvec)*1.2]);
set(gca,'ylim',[0,max(DIFF_cmpts)*1.1]);




