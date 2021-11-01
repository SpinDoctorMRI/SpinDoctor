function lap_eig = load_laplace_eig(path, mf, diffusivity)
%LOAD_LAPLACE_EIG Load laplace eigendecomposition.
%
%   path: string - path to the lap_eig file or the folder containing lap_eig
%   mf (optional): struct
%   diffusivity: double
%
%   lap_eig: struct


lap_eig = [];
if nargin == 1
    % path points to the lap_eig file
    mfile = matfile(path, "Writable", false);
    lap_eig = mfile.lap_eig;
    for ilapeig = 1:length(lap_eig)
        % Check integrity of saved eigenfunctions
        if lap_eig(ilapeig).md5 ~= DataHash(lap_eig(ilapeig).funcs)
            mfile.lap_eig = [];
            lap_eig = [];
            break;
        end
    end
elseif nargin == 3
    % path points to the folder containing lap_eig
    eigenpath = path;
    filename = sprintf( ...
        "%s/lap_eig_lengthscale%g_neigmax%g.mat", ...
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
            [ls, nmax] = sscanf(lapeig_list(ilist).name, 'lap_eig_lengthscale%f_neigmax%f.mat');
            if ls <= length_scale && nmax >= neig_max
                % reuse larger eigendecomposition
                name = fullfile(lapeig_list(ilist).folder, lapeig_list(ilist).name);
                lap_eig = load_laplace_eig(name);
                if ~isempty(lap_eig)
                    fprintf("Rompute Laplace eigendecomposition using saved result: %s", name);
                    eiglim = length2eig(mf.length_scale, diffusivity);
                    lap_eig = reset_lapeig(lap_eig, eiglim, mf.neig_max);
                    break;
                end
            end
        end
    end
end
