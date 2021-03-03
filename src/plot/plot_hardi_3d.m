function plot_hardi_3d(directions, data, fig_title)
%PLOT_HARDI_3D Plot data in directions.
%
%   directions: struct
%   data: double(nsequence[, namplitude], 1)
%   fig_title


% Check if data is real
if ~isreal(data)
    warning("Data is complex. Plotting real part of data.");
    data = real(data);
end

% Extract HARDI points
points = directions.points;

% Sizes
if length(size(data)) == 2
    data = shiftdim(data, -1);
end
[namplitude, nsequence, ndirection] = size(data);

% Create surface triangulation
DT = delaunayTriangulation(points');
[c_sphere, ~] = convexHull(DT);

for iseq = 1:nsequence

    % Create figure title label if more than one experiment
    if nsequence == 1
        experiment_str = "";
    else
        experiment_str = sprintf(", experiment %d of %d", iseq, nsequence);
    end

    % Create figure title label if more than one b-value
    if namplitude == 1
        bvalue_str = "";
    else
        tmp = sprintf(repmat("%d ", 1, namplitude-1) + "%d", 1:namplitude);
        bvalue_str = sprintf(", b-values [%s] of %d", tmp, namplitude);
    end

    % Make figure
    figure;
    hold on;

    for iamp = 1:namplitude
        % Scale unit sphere by data
        sig = shiftdim(data(iamp, iseq, :), 2);
        sigpoints = points' .* sig;

        % Make surface plot
        % h = trisurf(c_sphere, sigpoints(:, 1), sigpoints(:, 2), sigpoints(:, 3));
        h = trisurf(c_sphere, sigpoints(:, 1), sigpoints(:, 2), sigpoints(:, 3), sig);
        if namplitude == 1
            set(h, "facealpha", 0.6);
        else
            % Add transparency to reveal all the layers
            set(h, "facealpha", 0.1);
        end
        set(h, "facecolor", "interp");
        set(h, "edgecolor", "none");

        % Plot nodes
        plot3(sigpoints(:, 1), sigpoints(:, 2), sigpoints(:, 3), "k.", "markersize", 2);

        % Make quiver plot
        % a = zeros(ndirection, 1);
        % quiver3(a, a, a, sig(:, 1), sig(:, 2), sig(:, 3));
    end

    % Title
    title(sprintf("%s, %d directions%s%s", ...
        fig_title, ndirection, experiment_str, bvalue_str));

    % Space dimensions
    axis equal;

    % Colorbar
    colorbar;
    % caxis([0, 1]);

    % Labels
    xlabel("x");
    ylabel("y");
    zlabel("z");
    grid on;

    % Set view
    view(3);
    % view([2, 2, 1]);
end
