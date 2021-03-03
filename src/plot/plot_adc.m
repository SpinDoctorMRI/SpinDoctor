function plot_adc(adc, adc_allcmpts, title_str)
%PLOT_ADC Plot ADC.
%
%   adc
%   adc_allcmpts
%   title_str (BTPDE/HADC/STA)
%
%   figure for ADC of nsequence experiments of (BTPDE/HADC/STA)

ncompartment = size(adc, 1);
nsequence = size(adc, 2);

xticklabel = num2cell(1:ncompartment);
xticklabel{ncompartment + 1} = "Total";

figure;
sgtitle(title_str);
for iseq = 1:nsequence
    subplot(nsequence, 1, iseq);
    hold on;
    bar(1:ncompartment, adc(:, iseq), "b");
    bar(ncompartment+1, adc_allcmpts(iseq), "r");
    title(sprintf("Experiment %d", iseq));
    set(gca, "xtick", 1:ncompartment + 1);
    set(gca, "xticklabel", xticklabel);
    % set(gca, "ylim", [0, max(adc_allcmpts)]);
    % set(gca, "ytick", linspace(0, max(adc_allcmpts), 6));
    xlabel("Compartment");
    ylabel("ADC");
    grid on;
end
