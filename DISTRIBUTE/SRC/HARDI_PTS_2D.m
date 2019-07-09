function [points_gdir,graddir_index,negii] = HARDI_PTS_2D(ngdir_total)

% compute the diffusion-encoding directions uniformly distributed on the circle.

thetavec = linspace(0,pi,ngdir_total+1);

points_gdir(:,1) = cos(thetavec(1:end-1)');
points_gdir(:,2) = sin(thetavec(1:end-1)');
points_gdir(:,3) = zeros(size(thetavec(1:end-1)'));

graddir_index = 1:ngdir_total;
ndir = length(graddir_index);
clear negii;
for idir = 1:ndir
    negii{idir} = [];
end
