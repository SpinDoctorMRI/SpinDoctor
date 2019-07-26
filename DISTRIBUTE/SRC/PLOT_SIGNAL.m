function PLOT_SIGNAL(bvalues,MF_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_allcmpts,title_str)

% plot the dMRI signal
% 
% Input:
%     1. bvalues
%     2. MF_allcmpts
%     3. Sig_free
%     4. ADC_allcmpts_S0
%     5. ADC_allcmpts
%     6. title_str
%     
% Output:
%     1 figure for the signal with title of title_str

markervec_cell = cell(6,1);
markervec_cell{1} = 'o';
markervec_cell{2} = 'x';
markervec_cell{3} = '+';
markervec_cell{4} = 'd';
markervec_cell{5} = 's';
markervec_cell{6} = '*';

linevec_cell = cell(6,1);
linevec_cell{1} = '-.';
linevec_cell{2} = ':';
linevec_cell{3} = '-.';
linevec_cell{4} = '-';
linevec_cell{5} = '--';
linevec_cell{6} = ':';

nexperi = size(bvalues,1);

figure; hold on
iplot = 0;

for iexperi = 1:nexperi   
    yvec = real(MF_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        ['b',markervec_cell{mod(iexperi-1,6)+1}]);
    set(h,'MarkerSize', 9, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = ['Experi ',num2str(iexperi)];
end

for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);
    yvec = ADC_allcmpts_S0(iexperi)*exp(-ADC_allcmpts(iexperi)*bvec);
    h = plot(bvec, log10(yvec),...
        ['b',linevec_cell{mod(iexperi-1,6)+1}]);
    set(h,'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = ['ADC fit. Experi ',num2str(iexperi)];
end

yvec = Sig_free;
plot(bvalues(:), log10(yvec), ['r','-']);
iplot = iplot + 1;
legend_vec{iplot} = ['free diffusion'];

legend(legend_vec{1:iplot});
legend('Location','southwest');
set(gca,'FontSize',12);
xlabel('b-value (s/mm^2)');
ylabel('log10(Signal)')
grid on;
title(title_str);