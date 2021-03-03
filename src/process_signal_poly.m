function [adc, kur, S0] = process_signal_poly(data, bvalues, bmin, bmax)
%PROCESS_SIGNAL_POLY Fit polynomial to log-signal against the b-values.
%
% Inputs
%     1. data: double(1, namplitude)
%     2. bvalues: double(1, namplitude)
%     3. bmin: double
%     4. bmax: double
%
% Outputs
%     1. adc: double
%     2. kur: double
%     3. S0: double


% Only consider b-values within range [bmin, bmax]
ind_b = (bvalues >= bmin & bvalues <= bmax);
bvalues = bvalues(ind_b);
data = data(ind_b);
logdata = log(data);

% Sizes
namplitude = size(bvalues, 2);

% Keep increasing degree until ADC is stable
found = false;
adc_old = 1;
kur_old = 1;
degree = 0;

% Check for special cases
if namplitude == 1
    error("Cannot fit data from one point.");
elseif namplitude == 2
    warning("Fitting ADC from two points. To assure precision, make sure " ...
        + "that both b-values are sufficiently small.")
    coeffs = polyfit(bvalues, logdata, 1);
    adc = -coeffs(1);
end

while ~found && degree <= namplitude - 2

    % Not converged yet. Increase polynomial degree
    degree = degree + 1;

    % Fit polynomial coefficients
    coeffs = polyfit(bvalues, logdata, degree);
    adc = -coeffs(degree);
    if degree == 1
        kur = 0;
    else
        kur = coeffs(degree - 1) / adc^2 * 6;
    end

    % Polynomial error
    diff = polyval(coeffs, bvalues) - logdata;

    % Compute convergence criteria
    converged = max(abs(diff)) <= 1e-3 * max(abs(logdata));
    adc_stable = abs(adc-adc_old) <= 1e-6 || abs(adc-adc_old) <= 0.05 * abs(adc_old);
    kur_stable = kur_old <= 0.15 ...
        || abs(kur-kur_old) < 0.15 * abs(kur_old)...
        || abs(kur-kur_old) < 0.15;

    % Check for convergence
    if converged && adc_stable && kur_stable
        found = true;
    end

    % Update values from previous iteration
    adc_old = adc;
    kur_old = kur;
end

% Check if the loop was stopped by reaching the highest polynomial degree
if ~found
    warning("Kurtosis may not be accurate.");
end

% Estimate signal for zero b-value
S0 = exp(coeffs(end));
