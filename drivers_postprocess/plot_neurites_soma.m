function plot_neurites_soma(femesh_soma,femesh_neurites,highlight_neurites)
%%PLOT_NEURITES_SOMA plots the soma and neurites meshes together
%
%   femesh_soma
%   fmesh_neurites
%   highlight_cmpt: (optional) indices of neurite to highlight.

% fig = figure;
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

if nargin == 2
center = mean(points, 2);
center(3) = max(points(3, :)) + 1;
text(center(1), center(2), center(3), "Soma", "color", "r", "fontsize", 24);
end
nneurites = length(femesh_neurites);
if nargin == 2
    cmap = colormap("winter");
    stepsize = floor(size(cmap,1)/nneurites);
    colors = cmap((1:nneurites)*stepsize,:);
    for ib = 1:nneurites
        femesh= femesh_neurites{ib};
        surface = triangulation([femesh.facets{:}]',[femesh.points{:}]');
        h = trisurf(surface);
        set(h, "facealpha", 0.25);
        set(h, "facecolor", colors(ib,:));
        set(h, "LineStyle", "none");
        points = [femesh.points{:}];
        pmin = min([pmin,points], [], 2);
        pmax = max([pmax,points], [], 2);
    
        center = mean(points, 2);
        center(3) = max(points(3, :)) + 1;
        text(center(1), center(2), center(3), sprintf("Neurite %d",ib), "color", colors(ib,:), "fontsize", 24);
    end
else
    for ib = 1:nneurites
        femesh = femesh_neurites{ib};
        surface = triangulation([femesh.facets{:}]',[femesh.points{:}]');
        h = trisurf(surface);
        points = [femesh.points{:}];
        pmin = min([pmin,points], [], 2);
        pmax = max([pmax,points], [], 2);
        if ismember(ib,highlight_neurites)
            alpha = 0.5; color = 'g';           
            center = mean(points, 2);
            center(3) = max(points(3, :)) + 1;
            text(center(1), center(2), center(3), sprintf("Neurite %d",ib), "color", 'k', "fontsize", 24);
        else
            alpha = 0.25; color = 'b';
        end
        set(h, "facealpha",alpha);
        set(h, "facecolor",color);
        set(h, "LineStyle", "none");
       
    
        
    end
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
