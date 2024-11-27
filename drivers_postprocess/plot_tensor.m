function plot_tensor(I,title_str)
figure;
[X,Y,Z]  = sphere(10);
sz = size(X);
X = reshape(X,1,[]); Y = reshape(Y,1,[]); Z = reshape(Z,1,[]);

% V = I \[X;Y;Z];
V = I*[X;Y;Z];
X = reshape(V(1,:),sz); Y = reshape(V(2,:),sz); Z = reshape(V(3,:),sz);

h = surf(X,Y,Z); hold on;

set(h, "facealpha", 0.8);
set(h, "facecolor", "interp");
set(h, "edgecolor", "none");
title('Inertia tensor')
grid on;
axis equal;
view(3)
if nargin ==2
    title(title_str);
end
if  min(eig(I))/max(eig(I)) < 0.01
    msg = sprintf("Tensor has large varition in eigenvalues:\n%f\n%f\n%f",sort(eig(I)));
    warning(msg)
end