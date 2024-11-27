function savepath = add_mf_str_savepath(savepath,mf)
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
end