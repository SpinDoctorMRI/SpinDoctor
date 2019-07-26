function PLOT_SIGNAL_BTPDE_MF(experi_btpde,SIG_allcmpts,SIG_MF_allcmpts,SIG_MFGA_allcmpts)

% Plot BTPDE, MF and MFGA signal
markervec_cell = cell(6,1);
markervec_cell{1} = 'o';
markervec_cell{2} = 'o';
markervec_cell{3} = '+';
markervec_cell{4} = 'd';
markervec_cell{5} = 's';
markervec_cell{6} = '*';

colorvec_cell{1} = 'r';
colorvec_cell{2} = 'b';
colorvec_cell{3} = 'k';
colorvec_cell{4} = 'c';
colorvec_cell{5} = 'g';
colorvec_cell{6} = 'm';

bvalues = experi_btpde.bvalues;
nexperi = size(bvalues,1);

figure;
hold on

%% BTPDE
 iplot = 0;
for iexperi = 1:nexperi   
    yvec = real(SIG_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{mod(iexperi-1,6)+1},markervec_cell{mod(0,6)+1}]);
    
    set(h,'MarkerSize', 9, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = ['BTPDE Experi ',num2str(iexperi)];
end
legend(legend_vec{1:iplot});
legend('Location','southwest');


%% MF
markervec_cell{1} = '-';
markervec_cell{2} = '-';
markervec_cell{3} = '--';
markervec_cell{4} = 'd';
markervec_cell{5} = 's';
markervec_cell{6} = 'x';
for iexperi = 1:nexperi   
    yvec = real(SIG_MF_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{mod(iexperi-1,6)+1},markervec_cell{mod(0,6)+6}]);
    set(h,'MarkerSize', 9, 'LineWidth',1);
    
    iplot = iplot + 1;
    legend_vec{iplot} = ['MF Experi ',num2str(iexperi)];
end
legend(legend_vec{1:iplot});
legend('Location','southwest');

%% MFGA
markervec_cell{1} = '--';
markervec_cell{2} = '--';
markervec_cell{3} = '--';
markervec_cell{4} = 'd';
markervec_cell{5} = 's';
markervec_cell{6} = '*';
for iexperi = 1:nexperi   
    yvec = real(SIG_MFGA_allcmpts(iexperi,:));
    bvec = bvalues(iexperi,:);
    h = plot(bvec, log10(yvec),...
        [colorvec_cell{mod(iexperi-1,6)+1},markervec_cell{mod(iexperi-1,6)+1}]);
    set(h, 'LineWidth',1);
    iplot = iplot + 1;
    legend_vec{iplot} = ['MFGA Experi ',num2str(iexperi)];
end

%%
legend(legend_vec{1:iplot});
legend('Location','southwest');
set(gca, 'FontSize',12)
grid on;
title('BTPDE, MF and MFGA signals');