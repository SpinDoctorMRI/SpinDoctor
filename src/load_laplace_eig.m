function lap_eig = load_laplace_eig(eigenpath, mf, diffusivity)
%LOAD_LAPLACE_EIG Load laplace eigendecomposition.
%
%   eigenpath: string - path to the lap_eig file or the folder containing lap_eig
%   mf: struct
%   diffusivity: double - mean diffusivity
%
%   lap_eig: struct


lap_eig = [];
if nargin == 1
    % path points to the lap_eig file
    fprintf("Load Laplace eigenfunctions from %s.\n", eigenpath);
    mfile = matfile(eigenpath, "Writable", false);
    lap_eig = mfile.lap_eig;
elseif nargin == 3
    % eigenpath points to the folder containing lap_eig
    filename = sprintf( ...
        "%s/lap_eig_lengthscale%.4f_neigmax%g.mat", ...
        eigenpath, mf.length_scale, mf.neig_max ...
    );

    if isfile(filename)
        % saved lap_eig exists
        lap_eig = load_laplace_eig(filename);
    end

    if isfolder(eigenpath) && isempty(lap_eig)
        % saved lap_eig doesn't exist, check existence of any larger lap_eig
        lapeig_list = dir(sprintf("%s/lap_eig_lengthscale*_neigmax*.mat", eigenpath));
        for ilist = 1:length(lapeig_list)
            ls_nmax = sscanf(lapeig_list(ilist).name, 'lap_eig_lengthscale%f_neigmax%f.mat');
            % check if there is a smaller length scale
            if mf.length_scale >= ls_nmax(1)
                % reuse larger eigendecomposition
                name = fullfile(lapeig_list(ilist).folder, lapeig_list(ilist).name);
                lap_eig = load_laplace_eig(name);
                if ~isempty(lap_eig)
                    fprintf("Compute Laplace eigendecomposition using saved result: %s\n", name);
                    eiglim = length2eig(mf.length_scale, diffusivity);
                    % reset larger lapeig according to eiglim and mf.neig_max
                    lap_eig = reset_lapeig(lap_eig, eiglim, mf.neig_max);
                    return;
                end
            end
        end
        % check if there is a larger neigmax (length scale is prior to neigmax)
        for ilist = 1:length(lapeig_list)
            ls_nmax = sscanf(lapeig_list(ilist).name, 'lap_eig_lengthscale%f_neigmax%f.mat');
            if ~isinf(ls_nmax(2)) && ls_nmax(2) >= mf.neig_max
                % reuse larger eigendecomposition
                name = fullfile(lapeig_list(ilist).folder, lapeig_list(ilist).name);
                lap_eig = load_laplace_eig(name);
                if ~isempty(lap_eig)
                    fprintf("Compute Laplace eigendecomposition using saved result: %s\n", name);
                    eiglim = length2eig(mf.length_scale, diffusivity);
                    lap_eig = reset_lapeig(lap_eig, eiglim, mf.neig_max);
                    return;
                end
            end
        end
    end
end
