function [x,y,z]= make_cylinder(centre,centre2,radius,radius2,num_points)
    N = max(floor(sqrt(num_points)), 3 );
    rprofile = linspace(radius,radius2,N);
    direction = centre2 - centre;
    height = norm(direction);
    [X,Y,Z] = cylinder(rprofile,N);
    R  = make_rotation_matrix(direction/norm(direction));
    Z = height*Z;
    points = R*([X(:),Y(:),Z(:)]');
    x = reshape(points(1,:),size(X))+centre(1);
    y = reshape(points(2,:),size(Y))+centre(2);
    z = reshape(points(3,:),size(Z))+centre(3);
end
