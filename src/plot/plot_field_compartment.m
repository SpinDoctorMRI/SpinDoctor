function plot_field_compartment(femesh, field, icmpt, title_str, ifield, focus_icmpt)
%PLOT_FIELD_COMPARTMENT Plot a field defined on the finite element mesh.
%   The field may be a magnetization, an eigenfunction...
%
% Only the field plot for the icmpt-th compartment is created. 
%
%   femesh: struct
%   field: cell(1, ncompartment)
%       Field to plot. For a given compartment, field{icmpt} is of the
%       form double(nnode(icmpt), nfield), where nfield is the number of
%       fields, usually 1. For neig eigenvalues, nfield=neig.
%   icmpt: int
%       Index of the compartment.
%   title_str: string
%       Figure title.
%   ifield (optional): int 
%       Index of field to plot. If not provided, it is set to 1.
%   focus_icmpt (optional): bool
%       If true, plot domain limits focus on the icmpt-th compartment.


% Check if user has provided a specific field column index
if nargin == (nargin(@plot_field_compartment) - 2)
    % The field has only one column
    ifield = 1;
    focus_icmpt = true;
elseif nargin < nargin(@plot_field_compartment)
    focus_icmpt = true;
end

% Determine limits of domain
if focus_icmpt
    pmin = min([femesh.points{icmpt}], [], 2);
    pmax = max([femesh.points{icmpt}], [], 2);
    axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];
else
    pmin = min([femesh.points{:}], [], 2);
    pmax = max([femesh.points{:}], [], 2);
    axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];
end

figure;
facets = [femesh.facets{icmpt, :}];
points = femesh.points{icmpt};

h = trisurf(facets', points(1, :), points(2, :), points(3, :), real(field{icmpt}(:, ifield)));
set(h, "facecolor", "interp");
set(h, "facealpha", 0.6);
set(h, "EdgeColor", "none");

axis equal;
axis(axis_vec);
colorbar("eastoutside");
xlabel("x");
ylabel("y");
zlabel("z");
grid on;
view(3)
title(sprintf("%s, compartment: [%d]", title_str, icmpt));
