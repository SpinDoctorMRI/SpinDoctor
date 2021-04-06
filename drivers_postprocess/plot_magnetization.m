% Plot direction averaged magnetizations for BTPDE and MF (the latter is
% computed). This plotting script may only be run after running the script
% `driver_mf` or `driver_mf_save_load`.


%% Plot mesh and emphasize important compartments
cmpts = [1 7];%[22 20 24 23 27];
% cmpts_interesting = [8 13 17 29];
title_str = "Compartments with high direction averaged magnetization";
plot_emphasized_compartments(femesh, cmpts, title_str);
set(gca, "fontsize", 14);

% exportgraphics(gca, "output/emphasized_compartments.png");


%% Plot BTPDE and MF magnetization
disp("Plotting the magnetization solution using BTPDE solver and using eigenfunctions expansion");

plot_btpde = 1;
plot_mf = 1;
plot_diff = 1;

% Gradient directions to plot (for average: 0)
dir_inds = 1;
% dir_inds = [17 49];
% dir_inds = [33];
% dir_inds = [64];
% dir_inds = [17 33 49 64];

maxkappa = max(setup.pde.permeability);

for iseq = 1:nsequence
    seq = setup.gradient.sequences{iseq};
    for iamp = 1:namplitude
        qvalue = setup.gradient.qvalues(iamp, iseq);
        bvalue = setup.gradient.bvalues(iamp, iseq);
        for idir = dir_inds
            % BTDPE and MF Magnetizations
            if idir == 0
                % Direction averaged magnetizations
                solpde = btpde.magnetization_avg(:, iamp, iseq);
                soleig = mf.magnetization_avg(:, iamp, iseq);
                title_str = "direction averaged magnetization";
            else
                % Directionalized magnetizations
                solpde = btpde.magnetization(:, iamp, iseq, idir);
                soleig = mf.magnetization(:, iamp, iseq, idir);
                title_str = sprintf("magnetization, dir=[%.2f; %.2f; %.2f]", ...
                    setup.gradient.directions.points(:, idir));
            end

            % % Title string
            % title_str = sprintf("%s\n\\kappa=%g", title_str, maxkappa);
            title_str = sprintf("\\kappa=%g", maxkappa);
            % if seq.delta == seq.Delta
            %     title_str = sprintf("%s, \\delta=\\Delta=%g", title_str, seq.delta);
            % else
            %     title_str = sprintf("%s, \\delta=%g, \\Delta=%g", title_str, seq.delta, seq.Delta);
            % end
            % title_str = sprintf("%s, q=%g", title_str, qvalue);
            % title_str = sprintf("%s, b=%g", title_str, setup.gradient.bvalues(iamp, iseq));
            % title_str = sprintf("%s, q=%g, b=%g", title_str, qvalue, setup.gradient.bvalues(iamp, iseq));
            gamma = 2.67513 * 1e-04;
            title_str = sprintf("%s, g=%.3f, b=%g", title_str, qvalue / gamma, setup.gradient.bvalues(iamp, iseq));

            % Max absolute value of BTPDE on the entire domain (for rel error)
            maxabs = max(cellfun(@(x) max(abs(x)), solpde));

            % Relative error field (relative to allcmpts max abs)
            relerr = cellfun(@(x, y) abs(x - y) ./ maxabs, soleig, solpde, "UniformOutput", false);

            % Take real parts for plotting
            % solpde = cellfun(@real, solpde, "UniformOutput", false);
            soleig = cellfun(@real, soleig, "UniformOutput", false);

            % % Determine bounds
            % minbtpde = min(cellfun(@min, solpde));
            % maxbtpde = max(cellfun(@max, solpde));
            % minmf = min(cellfun(@min, soleig));
            % maxmf = max(cellfun(@max, soleig));
            % cmin = min(minbtpde, minmf);
            % cmax = max(maxbtpde, maxmf);

            cmin = 0;
            cmax = 1;

            % Plot BTPDE
            if plot_btpde
                plot_field_everywhere(femesh, solpde, "BTPDE " + title_str);
                % plot_field(femesh, solpde, setup.pde.compartments, title_str)
                % caxis([cmin, cmax]);
                set(gca, "fontsize", 14);
                if idir
                    % Plot gradient direction
                    plot_gdir_arrow(setup.gradient.directions.points(:, idir));
                end
            end

            % Plot MF
            if plot_mf
                title_str = "MF " + title_str;
                plot_field_everywhere(femesh, soleig, title_str);
                % plot_field(femesh, soleig, setup.pde.compartments, "MF " + title_str)
                % caxis([cmin, cmax]);
                % set(gca, "fontsize", 14);
                if idir
                    % Plot gradient direction
                    plot_gdir_arrow(setup.gradient.directions.points(:, idir));
                end
            end

            % Plot difference between MF and BTPDE
            if plot_diff
                plot_field_everywhere(femesh, relerr, "Relative error MF BTPDE " + title_str);
                % plot_field(femesh, relerr, setup.pde.compartments, "MF " + title_str])

                % caxis([0,0.01]);% 1]);
                set(gca, "fontsize", 14);
                if idir
                    % Plot gradient direction
                    plot_gdir_arrow(setup.gradient.directions.points(:, idir));
                end
            end
            % colorbar("off");
            % view(-90, 90);
            % view(2);

            % title(title_str, "fontsize", 16);
            % fname = sprintf("output/mag/k%g_q%g_b%g.png", maxkappa, qvalue, bvalue);
            fname = sprintf("output/mag/k%g_q%g_b%g_dir%g.png", maxkappa, qvalue, bvalue, idir);
%             exportgraphics(gca, fname, "resolution", 300);
        end
    end
end

% Clear temporary variables
clear solpde soleig relerr


%% Plot compartment signals over volumes
if true
    signal = mf.signal;
    signal_allcmpts = mf.signal_allcmpts;
    title_str = "MF";
else
    signal = btpde.signal;
    signal_allcmpts = btpde.signal_allcmpts;
    title_str = "BTPDE";
end

weights = ones(ncompartment, 1);
weights(end) = weights(end) / 3;

% signal_allcmpts = squeeze(sum(weights .* signal(1:end, :, :, :), 1));
signal_allcmpts = squeeze(sum(signal(1:end-1, :, :, :), 1));


idir = 0; % 0 for avg, idir for direction idir
if idir == 0
    % Direction averaged signal
    signal = mean(signal, 4);
    signal_allcmpts = mean(signal_allcmpts, 3);
    dir_str = "direction averaged";
else
    % Directionalized signal
    signal = signal(:, :, :, idir);
    signal_allcmpts = signal_allcmpts(:, :, idir);
    dir_str = sprintf("dir=[%.2f; %.2f; %.2f]", directions.points(:, idir));
end

% Normalize by volume
% signal = signal ./ volumes';
signal = real(signal);
signal(signal < 0) = 0;
% signal = abs(signal);

% Sort volumes and exclude ECS (last compartment)
% [~, inds] = sort(volumes(1:end-1));
[~, inds] = sort(volumes(1:end));
% colors = {"#0072BD" "#D95319" "#EDB120" "#7E2F8E" "#77AC30" "#4DBEEE" "#A2142F"};
colors = ["r" "b" "k" "#008800" "m" "c" "g" "y"];
markers = ["x" "o" "d" "s" ">" "v" "^"];
linestyles = ["-" "-." "--" ":"];

maxkappa = max(setup.pde.permeability);

figure;
hold on;

a = get(gcf, "position");
set(gcf, "position", 0.7 * a);

labels = strings(1, 0);
for iseq = 1:nsequence
    for iamp = 1:namplitude
        seq = setup.gradient.sequences{iseq};
        qvalue = setup.gradient.qvalues(iamp, iseq);
        bvalue = setup.gradient.bvalues(iamp, iseq);

        % Label
%         labels(end + 1) = sprintf("b=%g, %s", bvalue, seq);
        labels(end + 1) = sprintf("b=%g", bvalue);

        signal_allcmpts = real(signal_allcmpts);

        % Plot signal over volumes
        % h = plot(volumes(inds) / sum(volumes(inds)), signal(inds, iamp, iseq) / signal_allcmpts(iamp, iseq));
        % h = plot(volumes(inds) / sum(volumes(inds)), signal(inds, iamp, iseq));
        % h = plot(volumes(inds) / sum(volumes(inds)), signal(inds, iamp, iseq) / sum(volumes));
        % h = plot(volumes(inds) / sum(volumes(inds)), signal(inds, iamp, iseq) ./ volumes(inds)');
        h = plot(volumes(inds), signal(inds, iamp, iseq) ./ volumes(inds)');
        % h = plot(sqrt(volumes(inds)), signal(inds, iamp, iseq) ./ volumes(inds)');
        h.Color = colors(iseq);
        h.Marker = markers(iseq);
        % h.LineStyle = linestyles(iseq);

        % % Plot free diffusion
        % hh = plot(volumes(inds) / sum(volumes(inds)), exp(-setup.pde.diffusivity(inds).' * bvalue));
        % hh.Color = colors(iseq);
        % hh.Marker = markers(iseq);
        % hh.LineStyle = linestyles(iseq);

        % set(gca, "yscale", "log");
    end
end

legend(labels, "location", "east");
xlabel("Volume");
% ylabel("Signal");
grid on

set(gca, "fontsize", 14)

gamma = 2.67513 * 1e-04;
% subtitle_str = sprintf("\\kappa=%g, q=%g, %s", maxkappa, qvalue, dir_str);
% subtitle_str = sprintf("\\kappa=%g, q=%g", maxkappa, qvalue);
subtitle_str = sprintf("\\kappa=%g, g=%.3f", maxkappa, qvalue / gamma);

% title(title_str + " compartment signals", subtitle_str);
% title(title_str + " compartment signals " + subtitle_str, "fontsize", 18);
title(subtitle_str, "fontsize", 18)

% xlim([0 volumes(inds(end-1))]);
% xlim([volumes(inds(1)) 200]);
% xlim([0 200]);
ylim([0 1]);
% xticks(0:20:200);
% set(gca, "YScale", "log");

% set(gca, "position", [0 0 0.5 0.5]);
% exportgraphics(gca, sprintf("output/k%g_q%g.png", maxkappa, qvalue));


% clear signal colors markers labels maxkappa idir


%% Functions
function plot_gdir_arrow(dir)
% Plot gradient direction
color = "r";
p = [24 29 0];
dir = 10 * dir;
q = quiver3(p(1), p(2), p(3), dir(1), dir(2), dir(3), 0);
q.LineWidth = 1.5;
q.Color = color;
q.MaxHeadSize = 1;
text(p(1), p(2)-3, p(3), "g", "fontsize", 30, "color", color, ...
    "VerticalAlignment", "baseline");
end
