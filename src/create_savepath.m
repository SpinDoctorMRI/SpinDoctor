function savepath = create_savepath(setup, simulation_method, root_path)
%CREATE_SAVEPATH Generate the folder path for saving simulation results.
%
%   setup: struct
%   simulation_method: string
%   root_path: root folder path
%
%   savepath: path string

if nargin < nargin(@create_savepath)
    root_path = 'saved_simul';
end

% Create folder for saving results
[~, name, ext] = fileparts(setup.name);
name = name + replace(ext, '.', '_');

% Encode geometry information
refinement_str = "";
if isfield(setup.geometry, "refinement")
    refinement_str = sprintf("_refine%g", setup.geometry.refinement);
end
ecs_str = sprintf("_%s", setup.geometry.ecs_shape);
if isfield(setup.geometry, 'ecs_ratio') && setup.geometry.ecs_shape ~= "no_ecs"
    ecs_str = ecs_str + sprintf("%g", setup.geometry.ecs_ratio);
end

if ismember(simulation_method, ["btpde", "btpde_midpoint"])
    save_dir_path_geometry = fullfile(root_path, name + refinement_str + ecs_str);

    % Encode pde parameters
    % Make sure different settings are saved in different folders using hash
    fields = {'initial_signal', 'mean_diffusivity'};
    pde = rmfields(setup.pde, fields);
    pde_str = "kappa";
    if isfield(pde, 'permeability_in_out')
        pde_str = pde_str + "_inout" + num2str(pde.permeability_in_out);
    end
    if isfield(pde, 'permeability_out_ecs')
        pde_str = pde_str + "_outecs" + num2str(pde.permeability_out_ecs);
    end
    pde_str = pde_str + sprintf("_out%g_md5_%s", ...
        pde.permeability_out, DataHash(pde, 10));

    % Construct the final path
    savepath = fullfile(save_dir_path_geometry, pde_str, simulation_method);
elseif ismember(simulation_method, ["lap_eig", "mf"])
    save_dir_path_geometry = fullfile(root_path, name + refinement_str + ecs_str);

    % Encode pde parameters
    % Laplace eigendecomposition doesn't depend on relaxation and initial density
    fields = {'relaxation_in', 'relaxation_out', 'relaxation_ecs', ...
    'initial_density_in', 'initial_density_out', 'initial_density_ecs', ...
    'initial_signal', 'mean_diffusivity'};
    pde = rmfields(setup.pde, fields);
    pde_str = "kappa";
    if isfield(pde, 'permeability_in_out')
        pde_str = pde_str + "_inout" + num2str(pde.permeability_in_out);
    end
    if isfield(pde, 'permeability_out_ecs')
        pde_str = pde_str + "_outecs" + num2str(pde.permeability_out_ecs);
    end
    pde_str = pde_str + sprintf("_out%g_md5_%s", pde.permeability_out, DataHash(pde, 10));

    % construct the final path
    if simulation_method == "mf"
        fields = {'initial_signal', 'mean_diffusivity'};
        pde = rmfields(setup.pde, fields);
        if isfield(pde, 'relaxation_in') && ~isfield(pde, 'relaxation_ecs')
            mf_str = sprintf("mf_relax_out%g_in%g_density_out%g_in%g_md5_%s", ...
                pde.relaxation_out, pde.relaxation_in, pde.initial_density_out, ...
                pde.initial_density_in, DataHash(pde, 10));
        elseif isfield(pde, 'relaxation_ecs') && ~isfield(pde, 'relaxation_in')
            mf_str = sprintf("mf_relax_out%g_ecs%g_density_out%g_ecs%g_md5_%s", ...
                pde.relaxation_out, pde.relaxation_ecs, ...
                pde.initial_density_out, pde.initial_density_ecs, DataHash(pde, 10));
        elseif isfield(pde, 'relaxation_ecs') && isfield(pde, 'relaxation_in')
            mf_str = sprintf("mf_relax_i%g_o%g_e%g_density_i%g_o%g_e%g_md5_%s", pde.relaxation_in, ...
                pde.relaxation_out, pde.relaxation_ecs, pde.initial_density_in, ...
                pde.initial_density_out, pde.initial_density_ecs, DataHash(pde, 10));
        else
            mf_str = sprintf("mf_relax_out%g_density_out%g_md5_%s", pde.relaxation_out, ...
                pde.initial_density_out, DataHash(pde, 10));
        end
        savepath = fullfile(save_dir_path_geometry, pde_str, mf_str);
    else
        savepath = fullfile(save_dir_path_geometry, pde_str);
    end

elseif ismember(simulation_method, "hadc")
    save_dir_path_geometry = fullfile(root_path, name + refinement_str + ecs_str);

    % Encode pde parameters
    fields = {'permeability_in_out', 'permeability_out_ecs', ...
    'permeability_in', 'permeability_out', 'permeability_ecs', ...
    'initial_signal', 'mean_diffusivity'};
    pde = rmfields(setup.pde, fields);
    
    if isfield(pde, 'relaxation_in') && ~isfield(pde, 'relaxation_ecs')
        pde_str = sprintf("relax_out%g_in%g_density_out%g_in%g_md5_%s", ...
            pde.relaxation_out, pde.relaxation_in, pde.initial_density_out, ...
            pde.initial_density_in, DataHash(pde, 10));
    elseif isfield(pde, 'relaxation_ecs') && ~isfield(pde, 'relaxation_in')
        pde_str = sprintf("relax_out%g_ecs%g_density_out%g_ecs%g_md5_%s", ...
            pde.relaxation_out, pde.relaxation_ecs, ...
            pde.initial_density_out, pde.initial_density_ecs, DataHash(pde, 10));
    elseif isfield(pde, 'relaxation_ecs') && isfield(pde, 'relaxation_in')
        pde_str = sprintf("relax_i%g_o%g_e%g_density_i%g_o%g_e%g_md5_%s", pde.relaxation_in, ...
            pde.relaxation_out, pde.relaxation_ecs, pde.initial_density_in, ...
            pde.initial_density_out, pde.initial_density_ecs, DataHash(pde, 10));
    else
        pde_str = sprintf("relax_out%g_density_out%g_md5_%s", pde.relaxation_out, ...
            pde.initial_density_out, DataHash(pde, 10));
    end

    % construct the final path
    savepath = fullfile(save_dir_path_geometry, pde_str, simulation_method);

else
    savepath = fullfile(root_path, name, simulation_method);
end

if ~isfolder(savepath)
    mkdir(savepath);
end
end
