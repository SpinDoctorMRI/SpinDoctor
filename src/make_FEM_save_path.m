function FEM_save_path = make_FEM_save_path(pde,mf,multi_lap_eig,savepath)
    % Folder for saving
    mf_str = sprintf("neig%g_ls%.4f", ...
        mf.neig_max, mf.length_scale);
    if mf.surf_relaxation
        mf_str = "surf_relaxation_" + mf_str;
    end
    if mf.single
        mf_str = mf_str + "_single";
    end
    if ~isinf(mf.neig_max)
        % if neig_max is inf, mf.eigs doesn't exist or is removed.
        mf_str = mf_str + sprintf("_%s", DataHash(mf.eigs, 6));
    end
    savepath = fullfile(savepath, mf_str);
    
    no_relaxation = all(isinf(pde.relaxation));
    zero_permeability = all(pde.permeability==0);
    if multi_lap_eig
        msg = "Currently multi-compartment meshes are not supported";
        error(msg);
    end
    

    % if no_relaxation
    %     % Set paths for finding matrices.
    %      relaxation_str = "no_relaxation_";
    % else
    %     relaxation_str = sprintf("relaxation%f_",pde.relaxation);
    % end

    % if zero_permeability
    %     permeability_str = "zero_permeability_";
    % else
    %     permeability_str = sprintf("permeability%f_",pde.permeability);
    % end

    % Set paths for finding matrices.

    FEM_save_path = sprintf("%s/FEM_lap_eig.mat", ...
                savepath);
end