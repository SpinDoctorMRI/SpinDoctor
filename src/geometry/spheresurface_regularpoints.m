function [points, facets, volume] = spheresurface_regularpoints(r, npoint_request, use_exact)
%SPHERESURFACE_REGULARPOINTS Create npoint_request evenly distributed points on the sphere of radius r.


% Use exact number of points unless user has untoggled the option
if nargin < nargin(@spheresurface_regularpoints)
    use_exact = true;
end

% Handle case with 2 or fewer points
if npoint_request < 4
    if use_exact
        error("npoint_request must be at least 4.");
    else
        npoint_request = 4;
    end
end


% Keep reducing approximate number of points until the exact number of points is
% smaller or equal to the requested number of points
npoint_approx = npoint_request;
while true
    a = 4 * pi / npoint_approx;
    d = sqrt(a);
    Mnu = round(pi / d);
    dnu = pi / Mnu;
    dphi = a / dnu;

    % Create unit sphere of approximately npoint_approx points
    points = [];
    ncount = 1;
    for m = 0:Mnu-1
        nu = pi * (m + 0.5) / Mnu;
        Mphi = round(2 * pi * sin(nu) / dphi);
        n = 0:Mphi-1;
        phi = 2 * pi / Mphi * n';
        points(ncount + n, 1:3) = [sin(nu)*cos(phi), sin(nu)*sin(phi), cos(nu)*ones(Mphi, 1)];
        ncount = ncount + Mphi;
    end

    % Resulting number of points
    npoint = size(points, 1);

    % Break loop if fewer or equal to desired number of points
    if ~use_exact || npoint <= npoint_request
        break
    end

    % Otherwise: try with fewer points
    npoint_approx = npoint_approx - 1;
end

if use_exact && npoint < npoint_request
    % There are not enough points. Add (npoint_request-npoint) random points of
    % unit norm
    newpoints = rand(npoint_request - npoint, 3);
    points(npoint + 1:npoint_request, 1:3) = newpoints ./ vecnorm(newpoints, 2, 2);
end

% Scale unit sphere to given radius
points = r * points;

DT = delaunayTriangulation(points);
[facets, volume] = convexHull(DT);
