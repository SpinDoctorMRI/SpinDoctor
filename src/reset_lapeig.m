function lap_eig = reset_lapeig(lap_eig, eiglim, neigmax)
%RESET_LAPEIG Reset Laplace eigenvalues, functions based on a new eiglim or neigmax.
%   lap_eig: struct
%   eiglim: double
%   neigmax: int
%
%   lap_eig: struct with fields
%       values: [neig x 1]
%       funcs: [npoint x neig]
%       length_scale: [neig x 1]
%       totaltime: [1 x 1]


for ilapeig = 1:length(lap_eig)
    % Load eigenvalues
    values = lap_eig(ilapeig).values;

    % Remove eigenvalues above interval defined by length scale
    inds_keep = values <= eiglim;
    if ~isinf(neigmax)
        inds_keep(neigmax+1:end) = 0;
    end
    % Reset lap_eig
    lap_eig(ilapeig).values = values(inds_keep);
    lap_eig(ilapeig).funcs = lap_eig(ilapeig).funcs(:, inds_keep);
    lap_eig(ilapeig).length_scales = lap_eig(ilapeig).length_scales(inds_keep);
end

% Display function evaluation time
disp("Done with reset eigendecomposition.");
