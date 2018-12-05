xmin = 0;
xmax=0;
ymin=0;
ymax=0;
zmin=0;
zmax=0;
for ict = 1:Ncmpt
    xx=max(Pts_cmpt_reorder{icmpt}(1,:));
    yx=max(Pts_cmpt_reorder{icmpt}(2,:));
    zx=max(Pts_cmpt_reorder{icmpt}(3,:));
    xn=min(Pts_cmpt_reorder{icmpt}(1,:));
    yn=min(Pts_cmpt_reorder{icmpt}(2,:));
    zn=min(Pts_cmpt_reorder{icmpt}(3,:));
    xmin = min(xmin,xn);
    xmax = max(xmax,xx);
    ymin = min(ymin,yn);
    ymax = max(ymax,yx);
    zmin = min(zmin,zn);
    zmax = max(zmax,zx);
end

figure;
cmptvec = [1:min(5,Ncmpt-1),Ncmpt];
for ict = 1:length(cmptvec)
    icmpt = cmptvec(ict);
    subplot(2,3,mod(ict-1,length(cmptvec))+1);
    Fac = [];
    for iboundary = 1:Nboundary
        Fac = [Fac,Fac_boundary_reorder{icmpt}{iboundary}];
    end
    h = trisurf(Fac',Pts_cmpt_reorder{icmpt}(1,:),Pts_cmpt_reorder{icmpt}(2,:),...
        Pts_cmpt_reorder{icmpt}(3,:),real(YOUT{end}{end}{icmpt}(:,end)));
    set(h,'facealpha',0.9);
    axis equal;
    axis([xmin,xmax,ymin,ymax,zmin,zmax]); colorbar;
    view(3);
    title('Solution,last experi,last b-value');
end

figure; hold on
iplot = 0;
for iexperi = 1:nexperi   
    yvec = real(MF_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{1},markervec_cell{iexperi}]);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = [' Experi ',mynum2str(iexperi)];
end
yvec = Sig_free*sum(IC_cmpts.*VOL);
h = plot(bvalues(:), log10(yvec),...
    [colorvec_cell{2},'','-']);
set(h,'MarkerSize', 10, 'LineWidth',1);
iplot = iplot + 1;
legend_vec{iplot} = ['free diffusion'];
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
ylabel('log10(Sig)')
title(['UG = [',mynum2str(UG),']']);

figure;
for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+1); hold on;
    bar(1:Ncmpt,[ADC(:,iexperi)],'b');
    bar(Ncmpt+1,ADC_allcmpts(iexperi,1),'r');
    title('ADC BT');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end

for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+2); hold on;
    bar(1:Ncmpt,ADC_PDE_formulation(:,iexperi),'b');
    bar(Ncmpt+1,ADC_PDE_allcmpts(iexperi,1),'r');
    title('ADC DE');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end

for iexperi = 1:nexperi
    subplot(nexperi,3,(iexperi-1)*3+3); hold on;
    bar(1:Ncmpt,ADC_STA(:,iexperi),'b');
    bar(Ncmpt+1,ADC_STA_allcmpts(iexperi,1),'r');
    title('ADC STA');
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
    grid on;
end

figure; 
subplot(2,3,1);
spy(boundary_mat); 
xlabel('iboundary'); ylabel('icmpt');
title('Connections boundary-compartment');
set(gca,'Ytick',[1:Ncmpt]);
set(gca,'Xtick',[1:Nboundary]);
grid on;

subplot(2,3,2); hold on;
bar(1:Ncmpt,DIFF_cmpts,'b');
bar(Ncmpt+1,ADC_free_allcmpts,'r');
title('DIFF free');
set(gca,'ylim',[0,max(DIFF_cmpts)]);
set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt');
grid on;

subplot(2,3,3); hold on;
bar(1:Nboundary,kappa_vec,'b');
title('Permeability');
set(gca,'ylim',[0,max(kappa_vec)]);
set(gca,'Ytick',linspace(0,max(kappa_vec),6));
set(gca,'Xtick',[1:Nboundary]);
xlabel('iboundary');
grid on;

subplot(2,3,4); hold on;
bar(1:Ncmpt,VOL,'b');
bar(Ncmpt+1,VOL_allcmpts,'r');
title('VOL');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt');
grid on;

subplot(2,3,5); hold on;
bar(1:Ncmpt,SA,'b');
bar(Ncmpt+1,sum(SA),'r');
title('Surface Area');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt');
grid on;

subplot(2,3,6); hold on;
bar(1:Ncmpt,SAu,'b');
bar(Ncmpt+1,sum(SAu),'r');
title('SA in U_g');
set(gca,'Xtick',[1:Ncmpt+1]);
xlabel('icmpt');
grid on;