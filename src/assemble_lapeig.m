function [values, funcs] = assemble_lapeig(lap_eig, M, npoint_cmpts)
%ASSEMBLE_LAPEIG Assemble eigenvalues and eigenfunctions of all decoupled compartments.
%
%   lap_eig: indexed struct with fields
%       values: double(neig, 1)
%       funcs: double(npoint, neig)
%   M: mass matrix
%   npoints_cmpts: int(1, ncompartment).
%
%   values: [neig_allcmpts x 1]
%   funcs: [npoint_allcmpts x neig_allcmpts]

fprintf("Assemble eigendecomposition of %d compartments.\n", length(lap_eig));

npoint_allcmpts = sum(npoint_cmpts);
neig_cmpts = zeros(1, length(lap_eig));

% assemble eigenvalues
values = [];
for ilapeig = 1:length(lap_eig)
    values = [values; lap_eig(ilapeig).values];
    neig_cmpts(ilapeig) = length(lap_eig(ilapeig).values);
end

% assemble eigenfunctions
point_inds = cumsum([0 npoint_cmpts]);
get_point_inds = @(icmpt) point_inds(icmpt)+1:point_inds(icmpt+1);
eig_inds = cumsum([0 neig_cmpts]);
get_eig_inds = @(icmpt) eig_inds(icmpt)+1:eig_inds(icmpt+1);

funcs = zeros(npoint_allcmpts, length(values));
for ilapeig = 1:length(lap_eig)
    point_range = get_point_inds(ilapeig);
    eig_range = get_eig_inds(ilapeig);
    funcs(point_range, eig_range) = lap_eig(ilapeig).funcs;
end

% Sort eigenvalues in increasing order
[values, indices] = sort(values);
funcs = funcs(:, indices);

% Normalize eigenfunctions with mass weighting
funcs = funcs ./ sqrt(dot(funcs, M * funcs));
end
