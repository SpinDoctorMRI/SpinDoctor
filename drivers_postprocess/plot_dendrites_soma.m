function fig = plot_dendrites_soma(femesh_soma,dendrites)
%%PLOT_DENDRITES_SOMA plots the soma and dendrites meshes together
%
%   femesh_soma
%   dendrites
% 

fig = figure;
hold on;

warning('off','all');
soma = triangulation([femesh_soma.facets{:}]',[femesh_soma.points{:}]');
points = [femesh_soma.points{:}];
pmin = min(points, [], 2);
pmax = max(points, [], 2);

h = trisurf(soma);
set(h, "facealpha", 0.25);
set(h, "facecolor","r");
set(h, "LineStyle", "none");

center = mean(points, 2);
center(3) = max(points(3, :)) + 1;
text(center(1), center(2), center(3), "Soma", "color", "r", "fontsize", 24);

ndendrites = length(dendrites);
cmap = colormap("winter");
stepsize = floor(size(cmap,1)/ndendrites);
colors = cmap((1:ndendrites)*stepsize,:);

for iden = 1:ndendrites
    femesh_dendrite = dendrites{iden};
    dendrite = triangulation([femesh_dendrite.facets{:}]',[femesh_dendrite.points{:}]');
    h = trisurf(dendrite);
    set(h, "facealpha", 0.25);
    set(h, "facecolor", colors(iden,:));
    set(h, "LineStyle", "none");
    points = [femesh_dendrite.points{:}];
    pmin = min([pmin,points], [], 2);
    pmax = max([pmax,points], [], 2);

    center = mean(points, 2);
    center(3) = max(points(3, :)) + 1;
    text(center(1), center(2), center(3), sprintf("Dendrite %d",iden), "color", colors(iden,:), "fontsize", 24);

end
warning('on','all');

xlabel("x");
ylabel("y");
zlabel("z");
grid on;
axis equal;
axis_vec = [pmin(1) pmax(1) pmin(2) pmax(2) pmin(3) pmax(3)];
axis(axis_vec);
a = gca;
a.FontSize = 20;
view(3);
