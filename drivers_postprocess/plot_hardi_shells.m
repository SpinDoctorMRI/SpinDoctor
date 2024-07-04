function plot_hardi_shells(setup,signal)


directions= setup.gradient.directions;

data =squeeze(signal);data = real(data);
data = squeeze(data);
namplitude= setup.namplitude;

% Create surface triangulation
DT = delaunayTriangulation(directions');
[c_sphere, ~] = convexHull(DT);

sig_max = max(data,[],'all'); sig_min = min(data,[],'all');
for iamp = 1:namplitude
    % Scale unit sphere by data
    if namplitude > 1
        sig = data(iamp,:);
    else
        sig = data';
    end
    % 
    sigpoints = (directions.* sig)';

    % Make surface plot
    % h = trisurf(c_sphere, sigpoints(:, 1), sigpoints(:, 2), sigpoints(:, 3));
    h=trisurf(c_sphere, sigpoints(:, 1), sigpoints(:, 2), sigpoints(:, 3), sig);
    hold on;
    set(h, "facealpha", 0.5);
    set(h, "facecolor", "interp");
    set(h, "edgecolor", "none");
    % Plot nodes
    plot3(sigpoints(:, 1), sigpoints(:, 2), sigpoints(:, 3), "k.", "markersize", 5);

    
end
% Space dimensions
axis equal;
% clim([sig_min,sig_max]);
% Colorbar 
% caxis([0, 1]); 
% Labels
xlabel("x");
ylabel("y");
zlabel("z");
grid on;
view(3);
cb = colorbar(); 
title('Hardi plot')