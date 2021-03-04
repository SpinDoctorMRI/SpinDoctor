%% Save Laplace eigenvalues for given kappa (do this for all kappas)
idomain = 1;
lapeigvals_1eminusx = lap_eig.values{idomain};

save("output/lapeigvals_1eminusx.mat", "lapeigvals_1eminusx");

%% Plot Laplace eigenvalues for all kappas (assuming they are saved)
load output/lapeigvals_1eminus0;
load output/lapeigvals_1eminus2;
load output/lapeigvals_1eminus3;
load output/lapeigvals_1eminus4;
load output/lapeigvals_1eminus5;
load output/lapeigvals_1eminus9;

figure;
% a = get(gcf, "position");
% set(gcf, "position", 0.5*a);
hold on;

% style = {"linewidth", 1.5};
% style = {".", "marker", ".", "markersize", 10, "linewidth", 2};
style = {"-", "linewidth", 2};
% style = {};
plot(lapeigvals_1eminus0/mean_diffusivity, style{:});
plot(lapeigvals_1eminus2/mean_diffusivity, style{:});
plot(lapeigvals_1eminus3/mean_diffusivity, style{:});
plot(lapeigvals_1eminus4/mean_diffusivity, style{:});
plot(lapeigvals_1eminus5/mean_diffusivity, style{:});
plot(lapeigvals_1eminus9/mean_diffusivity, style{:});
maxval = max(lapeigvals_1eminus9) / mean_diffusivity;
maxeig = length(lapeigvals_1eminus9);
% lgd_cell = {"$\kappa=10^{-2}$", "$\kappa=10^{-3}$", "$\kappa=10^{-4}$", "$\kappa=10^{-5}$", "$\kappa=10^{-9}$"};
lgd_cell = {"$\kappa=1$", "$\kappa=10^{-2}$", "$\kappa=10^{-3}$", "$\kappa=10^{-4}$", "$\kappa=10^{-5}$", "$\kappa=10^{-9}$"};
lgd_cell{end + 1} = "$\lambda_{\max}$";
% line([2/3*maxeig maxeig+10], [maxval, maxval], "linestyle", "--", "color", "k");
line([2/3*maxeig maxeig+10], [maxval, maxval], "linewidth", 2, "linestyle", "--", "color", "k");
lgd = legend(lgd_cell);
set(lgd, "location", "northwest");
set(lgd, "interpreter", "latex");
set(lgd, "fontsize", 18);

xlabel("Eigenvalue index");
% ylabel("$\lambda/\sigma$", "Interpreter", "latex");
grid on
axis("tight");
% a = ylim;
% ylim([a(1) a(2)+0.03])
% ylim([0 0.17]);
% xlim([1 55]);
set(gca, "fontsize", 14);
title_str = "Laplace eigenvalues, L=3";
title(title_str, "fontsize", 18);

exportgraphics(gca, "output/lap.png");



%% Determine localization of BT eigenfunctions
eigindex_use_lap = 1:neig;
eigindex_use_btpde = 1:neig;
dirs_use = [1 20 40 55];
% support = zeros(namplitude, nsequence, ndirection, neig);
for iseq = 2%1:nsequence
    for iamp = 1%1:namplitude
        for idir = dirs_use
            clear bt_eig_funcs bt_eig_funcs_sep int_bt
            bt_eig_funcs(:, eigindex_use_btpde) = lap_eig.funcs{1}(:, eigindex_use_lap)...
                * bt_eig.invVsort{iamp, iseq, 1, idir}(eigindex_use_btpde, eigindex_use_lap)';

            figure;
            hold on;
            maxabs = 0;
            for ieig = eigindex_use_btpde
                flag = 0;
                bt_eig_funcs_sep = cell(1, ncompartment);
                for icmpt = 1:ncompartment
                    npoint_icmpt = size(femesh.points{icmpt}, 2);
                    bt_eig_funcs_sep{icmpt}(1:npoint_icmpt, :)...
                        = bt_eig_funcs(flag + 1:flag + npoint_icmpt, :);
                    flag = flag + npoint_icmpt;
                end

                % Mean value
                int_bt = integral_fem(bt_eig_funcs(:, ieig), femesh);
                int_bt = int_bt ./ volumes;

                % Determine max y-value
                maxabs = max(maxabs, max(abs(int_bt)));

                % Each eigenindex is assigned to the support compartment
                [~, support(iamp, iseq, idir, ieig)] = max(int_bt);

                % Sort in ascending order by volume
%                 [vol_sort, inds] = sort(volumes(1:end-1));
%                 plot(vol_sort, abs(int_bt(inds)), "-d");
            end
%             ylims = [-maxabs / 10, 1.1 * maxabs];
%             line([volumes(22) volumes(22)], ylims, "Color", [1 0 0]);
%             line([volumes(24) volumes(24)], ylims, "Color", [1 0 0]);
%
%             title("BT eigenfunctions, average value in each compartment, absolute value");
%             xlabel("Compartment volume, exluding ECS");
        end
    end
end


%% Plot complex BT eigenvalues
cmpts_boundary = [22 20 24 23 27];
cmpts_corner = [22 24];
cmpts_interesting = [8 13 17 29];

maxkappa = max(setup.pde.permeability);

markers = {"d", "+", "x", "o"};
imrk = 1;

figure;
hold on;
iseq = 2;
iamp = 1;

allvals = [];
for idomain = 1:ndomain
    for idir = dirs_use
        eigvals = bt_eig.Dsort{iamp, iseq, idir};
        eigvals = eigvals / mean_diffusivity;
        allvals = [allvals, eigvals(eigindex_use_btpde)];
        doplot = true;
        for ieig = eigindex_use_btpde
            eigval = eigvals(ieig);
            if support(iamp, iseq, idir, ieig) == ncompartment % ECS
                color = "g";
                % doplot = false;
            elseif any(support(iamp, iseq, idir, ieig) == cmpts_corner)
                color = "r";
            else
                color = "b";
            end
            if doplot
                plot(real(eigval), imag(eigval), "marker", markers{imrk}, "color", color);
                % text(real(eigval), imag(eigval), num2str(support(iamp, iseq, idir, ieig)), "color", color);
            end
        end
        imrk = imrk + 1;
    end
end
seq = setup.gradient.sequences{iseq};
qval = setup.gradient.qvalues(iamp, iseq);
experiment_str = sprintf("\\delta=%g, \\Delta=%g", seq.delta, seq.Delta);
bvalue_str = sprintf("q-value=%g", qval);
if length(dirs_use) == 1
    dir_str = sprintf("dir=[%.2f; %.2f; %.2f]", directions.points(:, dirs_use));
else
    dir_str = sprintf(["dirs [" + join(repmat("%d", 1, length(dirs_use))) + "]"], dirs_use);
end
title(sprintf("BT eigenvalues, \\kappa=%g, %s, %s, %s", maxkappa, experiment_str, bvalue_str, dir_str));
xlabel("Real part");
ylabel("Imaginary part");
grid on;

% valssort = sort(real(allvals));
% xlim([min(real(allvals)), valssort(400)])

% xlim([min(real(allvals)), max(real(allvals))])
% ylim([min(imag(allvals)), max(imag(allvals))])

%xlim = ([min(real]

%% Plot complex BT eigenvalues

figure;
hold on;
iseq = 1;
iamp = 1;
for idir = 1

    % The number of eigenfunctions
    neig = length(lap_eig.length_scales);

    eiguse = 1:neig;
    eigvals = bt_eig.Dsort{iamp, iseq, idir};

    eigvals = eigvals / diffusivity_domains;

    eigvals = eigvals(eiguse);

    plot(real(eigvals), imag(eigvals), "d");
end
title(sprintf("BT eigenvalues, gdir=[%.2f; %.2f; %.2f]", directions.points(:, idir)));
xlabel("Real part");
ylabel("Imaginary part");
