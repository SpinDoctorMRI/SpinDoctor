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

    % Remove redundant eigenvalues defined by maximum number of eigenvalues
    if ~isinf(neigmax)
        if ~isempty(inds_keep(neigmax+1:end)) && any(inds_keep(neigmax+1:end), 'all')
            msg = join([ ...
                "No eigenvalues were outside the interval.", ...
                "Consider increasing neig_max if there are more eigenvalues", ...
                "that may not have been found in the interval."
            ]);
            warning(msg);
        end
        inds_keep(neigmax+1:end) = false;
    end

    % Reset lap_eig
    lap_eig(ilapeig).values = values(inds_keep);
    lap_eig(ilapeig).funcs = lap_eig(ilapeig).funcs(:, inds_keep);
    lap_eig(ilapeig).length_scales = lap_eig(ilapeig).length_scales(inds_keep);
end
