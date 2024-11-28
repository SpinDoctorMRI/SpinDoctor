function plot_tensor(I,title_str)
figure;
if  min(eig(I))/max(eig(I)) < 0.01
    msg = sprintf("Tensor has large varition in eigenvalues:\n%f\n%f\n%f",sort(eig(I)));
    warning(msg)
    [V,E] = eig(I);
    [E,ind] = sort(diag(E));
    
    if E(2)/E(3) < 0.01 
        disp("Plotting as a stick");
        V = V(:,ind); V = V(:,3);
        t = linspace(-E(3),E(3),3);
        plot3(V(1)*t,V(2)*t,V(3)*t,LineWidth=5);
        shape = "stick";
    else 
        disp("Plotting as a disk")
        U = E(2)*V(:,2); V = E(3)*V(:,3);
        theta = linspace(0,2*pi,32);
        r = linspace(0,1,10);
        X = r'*cos(theta)*U(1) + r'*sin(theta)*V(1);
        Y = r'*cos(theta)*U(2) + r'*sin(theta)*V(2);
        Z = r'*cos(theta)*U(3) + r'*sin(theta)*V(3);
        h = surf(X,Y,Z); hold on;
        
        set(h, "facealpha", 0.8);
        set(h, "facecolor", "interp");
        set(h, "edgecolor", "none");
        axis equal;
        shape = "disk";
    end
else
    [X,Y,Z]  = sphere(10);
    sz = size(X);
    X = reshape(X,1,[]); Y = reshape(Y,1,[]); Z = reshape(Z,1,[]);
    
    V = I*[X;Y;Z];
    X = reshape(V(1,:),sz); Y = reshape(V(2,:),sz); Z = reshape(V(3,:),sz);
    
    h = surf(X,Y,Z); hold on;
    
    set(h, "facealpha", 0.8);
    set(h, "facecolor", "interp");
    set(h, "edgecolor", "none");
    axis equal;
    shape = "ellipsoid";
end
title(sprintf('Inertia tensor, shape = %s',shape))
grid on;
view(3)
if nargin ==2
    title(title_str);
end