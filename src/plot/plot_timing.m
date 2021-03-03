function plot_timing(elapsed_time, femesh, title_str, xlabel_str)
%PLOT_TIMING Plot elapsed time.
%
%   elapsed_time
%   femesh: struct
%
%   1 figure of computational times of nsequence experiments

npoint = sum(cellfun(@(x) size(x, 2), femesh.points));
nelement = sum(cellfun(@(x) size(x, 2), femesh.elements));

nsequence = size(elapsed_time, 2);
n1 = size(elapsed_time, 1);

xticklabel = num2cell(1:n1);
xticklabel{end + 1} = "Total";

figure;
sgtitle(title_str);
for iseq = 1:nsequence
    subplot(nsequence, 1, iseq);
    hold on;
    bar(1:n1, elapsed_time(:, iseq), "b");
    bar(n1+1, sum(elapsed_time(:, iseq)), "r");
    set(gca, "xtick", 1:n1+1);
    set(gca, "xticklabel", xticklabel);
    title(sprintf("Experi %d: %d nodes, %d elements", iseq, npoint, nelement))
    xlabel(sprintf("%s", xlabel_str));
    ylabel("Computational time (s)");
end
