function PLOT_TIMING(elapsed_time)

nexperi = size(elapsed_time,2);
nb = size(elapsed_time,1);
figure;
for iexperi = 1:nexperi
    subplot(nexperi,1,iexperi); hold on;
    bar(1:nb,[elapsed_time(:,iexperi)],'b');
	bar(nb+1,[sum(elapsed_time(:,iexperi))],'r');
    title(['Experi ',num2str(iexperi)]);
    ylabel('Computational time (s)');
    xlabel('Index (last: total)');
end