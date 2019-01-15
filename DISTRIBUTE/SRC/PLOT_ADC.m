function PLOT_ADC(ADC_cmpts,ADC_allcmpts,DIFF_cmpts,title_str)

nexperi = length(ADC_allcmpts);
Ncmpt = length(DIFF_cmpts);

figure;
for iexperi = 1:nexperi
    subplot(nexperi,1,iexperi); hold on;
    bar(1:Ncmpt,[ADC_cmpts(:,iexperi)],'b');
    bar(Ncmpt+1,ADC_allcmpts(iexperi,1),'r');
    title([title_str,' Experi ',num2str(iexperi)]);
    set(gca,'ylim',[0,max(DIFF_cmpts)]);
    set(gca,'Ytick',linspace(0,max(DIFF_cmpts),6));
	xlabel('icmpt (last: all cmpts)');
	ylabel('ADC_cmpts');
    
    grid on;
end



