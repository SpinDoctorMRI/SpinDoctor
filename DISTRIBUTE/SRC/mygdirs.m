function [points_gdir] = mygdirs(n)
% generate n directions evenly distributed in the x-y plan
    angle = linspace(0, 2*pi, n+1);
    angle = angle(1:end-1)';
    points_gdir = [cos(angle),sin(angle),zeros(n,1)];
end

