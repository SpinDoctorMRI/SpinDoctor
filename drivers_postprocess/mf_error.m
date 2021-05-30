%% Reference signal
tmp = split(save_dir_path_spindoctor, sprintf("refinement%g", setup.pde.refinement));
assert(length(tmp) == 2);
ref_str = tmp(1) + "refinement0.1" + tmp(end);
signal_ref = zeros(nsequence, namplitude, ndirection);
for iseq = 1:nsequence
    for iamp = 1:namplitude
        %%%% Extract experiment parameters
        seq = setup.gradient.sequences{iseq};
        bvalue = setup.gradient.bvalues(iamp, iseq);
        qvalue = setup.gradient.qvalues(iamp, iseq);

        %%%% Load data
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

%% MF signal
% [~, mf_signal_allcmpts, mf_itertimes] = compute_mf_signal_sub_bt(experiment, initial_signal_domains, lap_eig, bt_eig);
[~, mf_signal_allcmpts, mf_itertimes] = compute_mf_signal(experiment, initial_signal_domains, lap_eig);

%% Display computational time
disp("Laplace eigendecomposition:");
disp(lap_eig.itertimes);
% disp("BT eigendecomposition (per direction):");
% disp(bt_eig.itertimes / ndirection);
disp("MF computational time (per direction):");
disp(shiftdim(mean(mf_itertimes, 4), 1));

%%
disp("BTPDE computational time (per direction):");
disp(mean(btpde.itertimes, 3));

%% Error between MF and reference signal
meansig = mean(mf_signal_allcmpts, 3);
% meansig = mean(btpde.signal_allcmpts, 3);
% meansig_ref = mean(btpde.signal_allcmpts, 3);
meansig_ref = mean(signal_ref, 3);
err = abs(meansig - meansig_ref) ./ abs(meansig_ref);

% Display error in percent
disp("Error (%):");
fprintf("%.2f %.2f\n", 100 * err');


%%
% Number of BT eigenvalues
neig_bt = [1 10:10:300];


%% MF subset eig - signal
signal_allcmpts_reduce = cell(1, length(neig_bt));
for i = 1:length(neig_bt)
    [~, signal_allcmpts] = compute_mf_signal_sub_bt(experiment, initial_signal_domains, ...
        lap_eig, bt_eig, neig_bt(i));
    signal_allcmpts_reduce{i} = signal_allcmpts;
end


%% Plot relative deviation from full MF signal
colors = ["r" "b" "k" "#008800" "m" "c" "g" "y"];
markers = ["x" "o" "d" "s" ">" "v" "^"];
linestyles = ["-" "-." "--" ":"];
neig = length(lap_eig.values{1});
maxkappa = max(abs(setup.pde.permeability));

figure;

a = get(gcf, "position");
set(gcf, "position", 0.7 * a);

labels = strings(1, 0);

dir = 0;
if dir
    dir_str = sprintf("dir=[%.2f; %.2f; %.2f]", directions(:, dir));
else
    dir_str = "direction averaged";
end

for iseq = 1:nsequence
    seq = setup.gradient.sequences{iseq};
    for iamp  = 1:namplitude
        qvalue = setup.gradient.qvalues(iamp, iseq);
        bvalue = setup.gradient.bvalues(iamp, iseq);

        % Signals
        if dir
            % Directionalized
            signal_reduce = cellfun(@(x) x(iamp, iseq, dir), signal_allcmpts_reduce);
            signal_reference = btpde.signal_allcmpts(iamp, iseq, dir);
            % signal_reference = signal_ref(iamp, iseq, dir);
        else
            % Direction averaged
            signal_reduce = cellfun(@(x) mean(x(iamp, iseq, :)), signal_allcmpts_reduce);
            signal_reference = mean(btpde.signal_allcmpts(iamp, iseq, :));
            % signal_reference = mean(signal_ref(iamp, iseq, :));
        end

        % Error
        errorsig = abs(signal_reduce - signal_reference) / abs(signal_reference);

        % Plot
        h = semilogy(neig_bt, errorsig);
        set(h, "color", colors(iseq));
        set(h, "marker", markers(iseq));
        set(h, "linestyle", linestyles(iamp));
        hold on

        % Label
        if seq.delta == seq.Delta
            experiment_str = sprintf("\\delta=\\Delta=%g", seq.delta);
        else
            experiment_str = sprintf("\\delta=%g, \\Delta=%g", seq.delta, seq.Delta);
        end
        experiment_str = sprintf("b=%g", bvalue);
        labels(end+1) = experiment_str;
        % labels(end+1) = sprintf("%s, q=%g", experiment_str, qvalue);
    end
end
legend(labels);%, "location", "southwest");
xlabel("Number of eigenvalues");
% xlabel(sprintf("Number of Bloch-Torrey eigenvalues (%d Laplace eigenvalues)", neig));
% xlabel(sprintf("Number of Bloch-Torrey eigenvalues, sorted by real part_{ }", neig));
% xlabel(sprintf("Number of Bloch-Torrey eigenvalues, sorted by |V_{1j}|", neig));
% ylabel("Relative error");
set(gca, "fontsize", 14);
% axis("tight");
% xlim([0 200]);
ylim([1e-6 1]);
grid on

gamma = 2.67513 * 1e-04;
% title_str = "MF relative error, " + sprintf("\\kappa=%g, q=%g", maxkappa, qvalue)];
title_str = sprintf("\\kappa=%g, g=%.3f", maxkappa, qvalue / gamma);
title(title_str, "fontsize", 18);

exportgraphics(gca, sprintf("output/k%g_q%g.png", maxkappa, qvalue));
