function plot_field(femesh, field, compartments, title_str, ifield)
%PLOT_FIELD Plot a field defined on the finite element mesh.
%   The field may be a magnetization, an eigenfunction...
%
% Three figures are created, one for each domain of the mesh (IN, OUT, ECS). If
% one of these domains do not exist, the corresponding figure is not created.
%
%   femesh: struct
%   field: cell(1, ncompartment)
%       Field to plot. For a given compartment, field{icmpt} is of the
%       form double(nnode(icmpt), nfield), where nfield is the number of
%       fields, usually 1. For neig eigenvalues, nfield=neig.
%   compartments
%   title_str: Figure title.
%   ifield (optional): index of field to plot. If not provided, it is set to 1.


ncompartment = length(compartments);

% Check if user has provided a specific field column index
if nargin < nargin(@plot_field)
    % The field has only one column
    ifield = 1;
end

include_in = any(compartments == "in");
include_ecs = any(compartments == "ecs");

cmpts = "out";
figs.out = figure;
hold on
if include_in
    cmpts = ["in" cmpts];
    figs.in = figure;
    hold on
end
if include_ecs
    cmpts = [cmpts "ecs"];
    figs.ecs = figure;
    hold on
end


% Determine limits of domain
pmin = min([femesh.points{:}], [], 2);
pmax = max([femesh.points{:}], [], 2);
axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];

for icmpt = 1:ncompartment
    facets = [femesh.facets{icmpt, :}];
    points = femesh.points{icmpt};
    
    figure(figs.(compartments(icmpt)));
    h = trisurf(facets', points(1, :), points(2, :), points(3, :), real(field{icmpt}(:, ifield)));
    set(h, "facecolor", "interp");
    set(h, "facealpha", 0.6);
    set(h, "EdgeColor", "none");
end

for cmpt = cmpts
    figure(figs.(cmpt));
    axis equal;
    axis(axis_vec);
    colorbar("eastoutside");
    xlabel("x");
    ylabel("y");
    zlabel("z");
    grid on;
    view(3)
    inds = find(compartments == cmpt);
    cmpt_str = sprintf(join(repmat("%d", 1, length(inds))), inds);
    title(sprintf("%s, %s compartment: [%s]", title_str, cmpt, cmpt_str));
end
