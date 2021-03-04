function plot_signal(bvalues, signal_allcmpts, signal_free, adc_allcmpts_S0, adc_allcmpts, title_str)
%PLOT_SIGNAL Plot signal from solving the Bloch-Torrey equation.
%
%   bvalues
%      signal_allcmpts
%      signal_free
%      adc_allcmpts_S0
%      adc_allcmpts


markervec_cell = cell(6, 1);
markervec_cell{1} = "o";
markervec_cell{2} = "x";
markervec_cell{3} = "+";
markervec_cell{4} = "d";
markervec_cell{5} = "s";
markervec_cell{6} = "*";

nsequence = size(bvalues, 2);

figure;
hold on

iplot = 0;
for iseq = 1:nsequence
    yvec = real(signal_allcmpts(:, iseq)); % / real(signal_allcmpts(iseq, 1));
    bvec = bvalues(:, iseq);
    h = plot(bvec, log10(yvec), "b" + markervec_cell{mod(iseq - 1, 6) + 1});
    set(h, "MarkerSize",  10,  "LineWidth", 1);
    iplot = iplot + 1;
    legend_vec{iplot} = " Experi " + num2str(iseq);
end

yvec = signal_free(:);
h = plot(bvalues(:), log10(yvec), "r-");
set(h, "MarkerSize", 10, "LineWidth", 1);
iplot = iplot + 1;
legend_vec{iplot} = "free diffusion";
for iseq = 1:nsequence
    bvec = bvalues(:, iseq);
    yvec = adc_allcmpts_S0(iseq) * exp(-adc_allcmpts(iseq) * bvec);
    h = plot(bvec, log10(yvec), "b-");
    set(h, "MarkerSize", 10, "LineWidth", 1);
    iplot = iplot + 1;
    legend_vec{iplot} = sprintf("Experiment %d, linear part", iseq);
end

legend(legend_vec{1:iplot});
legend("Location", "southwest");
set(legend, "FontSize", 10)
set(gca, "FontSize", 10)
xlabel("B-value")
ylabel("log10(signal)")
grid on;
title(title_str);
