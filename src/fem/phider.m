function [dphi, detj, jac] = phider(coord, point, etype)
%PHIDER Return gradients of basis functions w.r.t. local coordinates (x, y, ...).
%   Copyright (c) 2013, Talal Rahman, Jan Valdman
%
% coord : coord(nod, nos, noe), the local coordinates of the
%         nodes which the shape functions are associated with.
% point : point(nod, nop), the coordinates of the
%         points on the reference element.
%  dphi : dphi(nod, nos, nop, noe), the gradients of the
%         shape functions (second) at all points (third)
%         with respect to the local cordinates.
%   jac : jac(nod, nod, nop, noe), the Jacobian matrices
%         at all nop points.
%  detj : detj(1, nop, noe), determinants of the Jacobian matrices
%         at all nop points
% etype : "P0", "P1", "P2", etc., the element type.
%
%         Note:
%         nod - dimension of the element.
%         nop - number of points.
%         nos - number of shape functions.
%         noe - number of elements.

jacout = nargout > 2;
detout = nargout > 1;

nod = size(coord, 1);
nop = size(point, 2);
nos = size(coord, 2);
noe = size(coord, 3);

% Derivatives with respect to the reference coordinates (xi, eta, ...)
dshape = shapeder(point, etype);

if jacout
    jac = zeros(nod, nod, nop, noe);
end
if detout
    detj = zeros(1, nop, noe);
end
dphi = zeros(nod, nos, nop, noe);

for poi = 1:nop
    tjac = smamt(dshape(:, :, poi), coord);
    [tjacinv, tjacdet] = aminv(tjac);
    dphi(:, :, poi, :) = amsm(tjacinv, dshape(:, :, poi));

    if jacout
        jac(:, :, poi, :) = tjac;
    end
    if detout
        detj(1, poi, :) = abs(tjacdet);
    end
end
