%% Reference signal
tmp = split(save_dir_path_spindoctor, sprintf("refinement%g", setup.pde.refinement));
assert(length(tmp) == 2);
ref_str = tmp(1) + "refinement0.1" + tmp(end);
signal_ref = zeros(nsequence, namplitude, ndirection);
for iseq = 1:nsequence
    for iamp = 1:namplitude
        % Extract experiment parameters
        bvalue = setup.gradient.bvalues(iamp, iseq);
        qvalue = setup.gradient.qvalues(iamp, iseq);
        seq = setup.gradient.sequences{iseq};

        % Load data
        if setup.gradient.values_type == "q"
            bvalue_str = sprintf("q%g", qvalue);
        else
            bvalue_str = sprintf("b%g", bvalue);
        end
        fname = sprintf("btpde_%s_%s_abstol%g_reltol%g.mat", class(seq), ...
            bvalue_str, setup.btpde.abstol, setup.btpde.reltol);
        call_cmd = sprintf("load %s/%s signal_allcmpts", ref_str, fname);
        disp(call_cmd);
        eval(call_cmd);

        % Store data
        signal_ref(iamp, iseq, :) = signal_allcmpts;
    end
end

% Display reference signal over zero signal
disp("Reference signal (S/S_0)");
disp(mean(signal_ref, 3) / sum(volumes));


%%
signal_allcmpts_reduce = btpde.signal_allcmpts;
maxkappa = max(abs(setup.pde.permeability));
fname = sprintf("output/btpde_sig_k%g_h%g.mat", maxkappa, setup.pde.refinement);
save(fname, "signal_allcmpts_reduce");

%% Number of Laplace eigenvalues
length_scales = 1:0.5:5;
n_length_scale = length(length_scales);
neig_lap = sum(lap_eig.length_scales{1} > length_scales);

%% MF subset eig - signal
% signal_allcmpts_reduce = cell(1, length(neig_lap));
for i = [2 4 6 8]%1:length(neig_lap)
    [~, signal_allcmpts] = compute_mf_signal(experiment, initial_signal, ...
        lap_eig, neig_lap(i));
    signal_allcmpts_reduce{i} = signal_allcmpts;
end

maxkappa = max(abs(setup.pde.permeability));
fname = sprintf("output/mf_sig_k%g_h%g.mat", maxkappa, setup.pde.refinement);
save(fname, "signal_allcmpts_reduce");

%%
maxkappa = max(abs(setup.pde.permeability));

load(sprintf("output/mf_sig_k%g_h0.5.mat", maxkappa));
signal_mf_h0p5 = signal_allcmpts_reduce;
load(sprintf("output/mf_sig_k%g_h0.2.mat", maxkappa));
signal_mf_h0p2 = signal_allcmpts_reduce;

load(sprintf("output/btpde_sig_k%g_h0.5.mat", maxkappa));
signal_btpde_h0p5 = signal_allcmpts_reduce;
load(sprintf("output/btpde_sig_k%g_h0.2.mat", maxkappa));
signal_btpde_h0p2 = signal_allcmpts_reduce;


%% Plot relative deviation from full MF signal
colors = ["r" "b" "k" "#008800" "m" "c" "g" "y"];
markers = ["x" "o" "d" "s" ">" "v" "^"];
linestyles = ["-" "-." "--" ":"];
maxkappa = max(abs(setup.pde.permeability));

dir = 0;
if dir
    % Directionalized
    dir_str = sprintf("dir=[%.2f; %.2f; %.2f]", directions(:, dir));
    signal_mf_h0p5 = cellfun(@(x) x(:, :, dir), signal_mf_h0p5, "UniformOutput", false);
    signal_mf_h0p2 = cellfun(@(x) x(:, :, dir), signal_mf_h0p5, "UniformOutput", false);
    signal_btpde_h0p5 = signal_btpde_h0p5(:, :, dir, :);
    signal_btpde_h0p2 = signal_btpde_h0p2(:, :, dir, :);
    signal_reference = signal_ref(:, :, dir);
else
    % Direction averaged
    dir_str = "direction averaged";
    signal_mf_h0p5 = cellfun(@(x) mean(x, 3), signal_mf_h0p5, "UniformOutput", false);
    signal_mf_h0p2 = cellfun(@(x) mean(x, 3), signal_mf_h0p2, "UniformOutput", false);
    signal_btpde_h0p5 = mean(signal_btpde_h0p5, 3);
    signal_btpde_h0p2 = mean(signal_btpde_h0p2, 3);
    signal_reference = mean(signal_ref, 3);
end

% Error
error_mf_h0p5 = abs(cat(3, signal_mf_h0p5{:}) - signal_reference) ./ abs(signal_reference);
error_mf_h0p2 = abs(cat(3, signal_mf_h0p2{:}) - signal_reference) ./ abs(signal_reference);
error_btpde_h0p5 = abs(signal_btpde_h0p5 - signal_reference) ./ abs(signal_reference);
error_btpde_h0p2 = abs(signal_btpde_h0p2 - signal_reference) ./ abs(signal_reference);

for iseq = 1:nsequence
    seq = setup.gradient.sequences{iseq};
    for iamp  = 1:namplitude
        qvalue = setup.gradient.qvalues(iamp, iseq);
        bvalue = setup.gradient.bvalues(iamp, iseq);

        figure;

        a = get(gcf, "position");
        set(gcf, "position", 0.7*a);

        % Plot

        semilogy(length_scales, shiftdim(error_mf_h0p5(iseq, iamp, :)), "color", "r", "marker", "o");
        hold on
        line([1 5], error_btpde_h0p5(iseq, iamp) * [1 1], "color", "r", "linestyle", "--");
        hold on

        semilogy(length_scales, shiftdim(error_mf_h0p2(iseq, iamp, :)), "color", "b", "marker", "d");
        hold on
        line([1 5], error_btpde_h0p2(iseq, iamp) * [1 1], "color", "b", "linestyle", "--");
        hold on

        % Label
        experiment_str = string(seq);

        legend("H=0.5, MF", "H=0.5, BTPDE", "H=0.2, MF", "H=0.2, BTPDE", "location", "northwest");
        xlabel("Minimum length scale");
        % ylabel("Relative error");
        set(gca, "fontsize", 14);
        % axis("tight");
        % xlim([0 200]);
        if iamp == 1
            ylim([1e-6 1]);
        else
            ylim([1e-4 1]);
        end
        grid on

        gamma = 2.67513 * 1e-04;
        % title_str = "MF relative error, " + sprintf("\\kappa=%g, q=%g" + maxkappa, qvalue);
        title_str = sprintf("\\kappa=%g, g=%.3f, \\delta=%g", maxkappa, qvalue / gamma, seq.delta);
        title(title_str, "fontsize", 18);

        fname = sprintf("output/mf_error/k%g_q%g_d%g.png", maxkappa, qvalue, seq.delta);
        exportgraphics(gca, fname, "resolution", 300);
    end
end
