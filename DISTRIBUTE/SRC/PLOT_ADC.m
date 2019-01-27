function PLOT_ADC(ADC_cmpts,ADC_allcmpts,DIFF_cmpts,title_str)

% plot ADC
% 
% Input:
%     1. ADC_cmpts
%     2. ADC_allcmpts
%     3. DIFF_cmpts
%     4. title_str (BTPDE/HADC/STA)
%     
% Output:
%     1 figure for ADC of nexperi experiments of (BTPDE/HADC/STA)

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
	ylabel('ADC\_cmpts');
    grid on;
end



