function plot_field_everywhere(femesh, field, title_str, ifield)
%PLOT_FIELD_EVERYWHERE Plot a field on all cmpts in the same figure.
%   The field may be a magnetization, an eigenfunction...
%
%   femesh: struct
%   field: cell(1, ncompartment)
%       Field to plot. For a given compartment, field{icmpt} is of the
%       form double(nnode(icmpt), nfield), where nfield is the number of
%       fields, usually 1. For neig eigenvalues, nfield=neig.
%   title_str: Figure title.
%   ifield (optional): index of field to plot. If not provided, it is set to 1.


% Check if user has provided a specific field column index
if nargin < nargin(@plot_field_everywhere)
    % The field has only one column
    ifield = 1;
end

ncompartment = femesh.ncompartment;

%% Determine limits of domain
pmin = min([femesh.points{:}], [], 2);
pmax = max([femesh.points{:}], [], 2);
axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];

%% Create figure
figure;
hold on
for icmpt = 1:ncompartment
    facets = [femesh.facets{icmpt, :}];
    points = femesh.points{icmpt};
    h = trisurf(facets', points(1, :), points(2, :), points(3, :), real(field{icmpt}(:, ifield)));
    set(h, "facealpha", 0.6);
    set(h, "facecolor", "interp");
    set(h, "edgecolor", "none");

    % % Plot nodes
    % plot3(points(1, :), points(2, :), points(3, :), ".");
end

%% Title
title(title_str);%, sprintf("All %d compartments", ncompartment));
% title(title_str, sprintf("All %d compartments", ncompartment));

%% Axis
    axis equal;
    axis(axis_vec);

%% Colorbar
colorbar("eastoutside");

%% Labels and grid
xlabel("x");
ylabel("y");
% zlabel("z");
% grid on;

%% View
view(3)
