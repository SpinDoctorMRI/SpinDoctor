function plot_hardi_2d(directions, signal, title_str)
%PLOT_HARDI_2D Plot direction dependent signal in the plane, for all experiments and b-values.
%
%   directions: struct
%   signal: double(nsequence, namplitude, ndirection) or double(nsequence, ndirection)


use_polar = true;

% Check if signal is real
if ~isreal(signal)
    warning("Signal is complex. Plotting real part of signal.");
    signal = real(signal);
end

% Check that the directions lie in the plane
if any(directions.points(3, :))
    warning("z components of gradient directions are not all zero." ...
        + " Projecting directions into the x-y plane.");
end

% Extract x and y components of gradient directions
points = directions.points(1:2, :);

% Sizes
if length(size(signal)) == 2
    signal = shiftdim(signal, -1);
end
[namplitude, nsequence, ndirection] = size(signal);

% Check if signal is also indexed for b-values
if namplitude == 1
    make_label = @(iamp, iseq) sprintf("Experiment %d", iseq);
else
    make_label = @(iamp, iseq) sprintf("Experiment %d, b-value %d", iseq, iamp);
end

% Angles
theta = cart2pol(points(1, :)', points(2, :)');

iplot = 1;
figure;

% Plot deformed unit circles
legend_vec = cell(1, namplitude * nsequence);
for iseq = 1:nsequence
    for iamp = 1:namplitude
        if use_polar
            sig = shiftdim(signal(iamp, iseq, :));
            h = polarplot(theta([1:end 1]), sig([1:end 1]));
            hold on
        else
            vec = zeros(ndirection + 1, 2);
            % Extract real part of signal
            sig = shiftdim(signal(iseq, iamp, :), 2);

            % Scale unit circle by signal
            vec(1:end-1, :) = sig .* points';

            % Close the deformed circle
            vec(end, :) = vec(1, :);

            % Plot signal scaled unit circle
            h = plot(vec(:, 1), vec(:, 2));
            hold on
        end
        set(h, "linewidth", 1.5);

        % Update legend
        legend_vec{iplot} = make_label(iamp, iseq);
        iplot = iplot + 1;
    end
end

% Plot origin
if ~use_polar
    origin = plot(0, 0, "x");
    axis equal;
    set(origin, "markersize", 12);
    set(origin, "linewidth", 1.5);
end

% Information
legend(legend_vec, "Location", "northwest");
grid on;
title(sprintf("%s in %d directions", title_str, ndirection));
