function PLOT_SIGNAL(bvalues,MF_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_allcmpts)

% plot signal from solving the Bloch-Torrey equation
% 
% Input:
%     1. bvalues
%     2. MF_allcmpts
%     3. Sig_free
%     4. ADC_allcmpts_S0
%     5. ADC_allcmpts
%     
% Output:
%     1 figure for the signal

markervec_cell = cell(6,1);
markervec_cell{1} = 'o';
markervec_cell{2} = 'x';
markervec_cell{3} = '+';
markervec_cell{4} = 'd';
markervec_cell{5} = 's';
markervec_cell{6} = '*';

nexperi = size(bvalues,1);

figure;
hold on
iplot = 0;
for iexperi = 1:nexperi   
    yvec = real(MF_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        ['b',markervec_cell{mod(iexperi-1,6)+1}]);
    set(h,'MarkerSize', 10, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = [' Experi ',num2str(iexperi)];
end
yvec = Sig_free;
h = plot(bvalues(:), log10(yvec),...
    ['r','','-']);
set(h,'MarkerSize', 10, 'LineWidth',1);
iplot = iplot + 1;
legend_vec{iplot} = ['free diffusion'];
for iexperi = 1:nexperi
    bvec = bvalues(iexperi,:);
    yvec = ADC_allcmpts_S0(iexperi)*exp(-ADC_allcmpts(iexperi)*bvec);
    h = plot(bvec, log10(yvec),...
        ['b','-']);
    set(h,'MarkerSize', 10, 'LineWidth',1);
end
legend(legend_vec{1:iplot});
set(legend, 'FontSize',10)
set(gca, 'FontSize',10)
xlabel('bvalue')
ylabel('log10(Sig)')