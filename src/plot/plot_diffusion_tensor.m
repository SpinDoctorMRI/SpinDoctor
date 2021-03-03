function plot_diffusion_tensor(diffusion_tensor, diffusivity)
%PLOT_DIFFUSION_TENSOR Plot diffusion tensor.
%
%   diffusion_tensor: double(3, 3, nsequence)
%   diffusivity: double


% Sizes
nsequence = size(diffusion_tensor, 3);

% Normalize diffusion tensor by diffusivity
diffusion_tensor = diffusion_tensor / diffusivity;

ymin = min(diffusion_tensor, [], "all");
ymax = max(diffusion_tensor, [], "all");

figure;
hold on;
sgtitle("Diffusion tensor")
for iseq = 1:nsequence
    subplot(1, nsequence, iseq);
    hold on;
    D = diffusion_tensor(:, :, iseq);
    yc = zeros(3, 1);
    yo = zeros(3, 1);
    yc(1:3) = D([1, 5, 9]);
    yo(1:3) = D([2, 3, 6]);
    bar(1:3, yc, "r");
    bar(4:6, yo, "b");
    % set(gca, "ylim", [min(0, min(dtensor, [], "all")), 1]);
    set(gca, "ylim", [ymin, ymax]);
    % set(gca, "ytick", linspace(0, 1, 3));
    set(gca, "xtick", 1:6);
    set(gca, "xticklabel", ["xx", "yy", "zz", "xy", "xz", "yz"]);
    % ylabel("");
    xlabel("Matrix index");
    title(sprintf("Experiment %d", iseq));
    grid on;
end
