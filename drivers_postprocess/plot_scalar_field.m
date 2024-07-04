function plot_scalar_field(v,s)
    
    DT = delaunayTriangulation(v');
    [f, ~] = convexHull(DT);
    s = reshape(s,[1,size(v,2)]);
    r =abs(s);
    h=trisurf(f, r.*v(1,:), r.*v(2,:), r.*v(3,:), real(s));
    set(h, "facealpha", 0.8);
    set(h, "facecolor", "interp");
    view(3);
    axis equal;
    % colorbar;
end