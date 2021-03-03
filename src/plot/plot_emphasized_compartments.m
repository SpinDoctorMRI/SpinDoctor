function plot_emphasized_compartments(femesh, cmpts, title_str)
%PLOT_EMPHASIZED_COMPARTMENTS Plot a field in all cmpts in the same figure.
%   The field may be a magnetization, an eigenfunction...
%
%   femesh: struct
%   cmpts: int(1, n_emphasized_compartments)
%   title_str: Figure title.


% Number of compartments
ncompartment = femesh.ncompartment;

%% Determine limits of domain
pmin = min([femesh.points{:}], [], 2);
pmax = max([femesh.points{:}], [], 2);
axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];

%% Create figure
figure;
hold on;
for icmpt = [ncompartment 1:ncompartment-1] % plot ECS first
    facets = [femesh.facets{icmpt, :}];
    points = femesh.points{icmpt};
    h = trisurf(facets', points(1, :), points(2, :), points(3, :));
    % set(h, "facealpha", 0.6);
    set(h, "facecolor", "interp");
    if ismember(icmpt, cmpts)
        set(h, "facecolor", "r");
        % set(h, "edgealpha", 0.2);
        % set(h, "edgecolor", "interp");
    elseif icmpt == ncompartment
        % ECS
        set(h, "edgealpha", 0.1);
        % set(h, "facealpha", 0.9);
        % set(h, "edgecolor", "none");
    else
        % set(h, "facealpha", 0.6);
        % set(h, "edgecolor", "interp");
        set(h, "edgealpha", 0.2);
        % set(h, "facecolor", "y");
    end

    % Plot nodes
    % plot3(points(1, :), points(2, :), points(3, :), ".");
end

%% Title
%title(sprintf("%s all compartments", title_str));
title(title_str);

%% Axis
axis equal;
axis(axis_vec);

%% Labels and grid
xlabel("x");
ylabel("y");
% zlabel("z");
grid on;

%% Set view
view(-65, 50); % Figures in paper
% view(2);
% view(3);
% view(-40, 60);
% view(90, 90);
