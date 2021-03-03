
zoom = false;
plot_femesh_everywhere(femesh, "");
x = [-12 7];
y = [2 8];
if zoom
    xlim(x);
    ylim(y);
    title(sprintf("H = %g", params_domain.refinement), "fontsize", 18);
else
    line([x(1) x(2) x(2) x(1) x(1)], [y(1) y(1) y(2) y(2) y(1)], 0.5*ones(1,5), "color", "r", "linewidth", 2);
    title(sprintf("Finite element mesh, H = %g", params_domain.refinement), "fontsize", 18);
end

exportgraphics(gca, sprintf("output/%g.png", params_domain.refinement));
