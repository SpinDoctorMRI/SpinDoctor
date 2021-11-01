%DRIVER_SAVE_LOAD Solve BPTDE, HADC and MF and save or load results.
%
%   It is highly recommended to read this driver to understand the workflow
%   of SpinDoctor.
%
%   The user is advised to read the latest version
%   from \url{https://github.com/jingrebeccali/SpinDoctor}

clear
restoredefaultpath

% Add SpinDoctor
addpath(genpath("src"));


%% Define inputs
% Get setup
addpath setups

% setup_1axon_analytical;
% setup_1sphere_analytical;
% setup_15spheres;
% setup_2axons_deform;
% setup_5axons_myelin_relax;
% setup_neuron;
% setup_4axons_flat;
% setup_30axons_flat;
setup_30axons;
% setup_200axons;

% Choose whether to save magnetization field
magnetization_flag = true;

% Prepare simulation
[setup, femesh, surfaces]  = prepare_simulation(setup);


%% Solve BTPDE
if isfield(setup, "btpde")
    disp("Computing or loading the BTPDE signals");
    savepath = create_savepath(setup, "btpde");
    btpde = solve_btpde(femesh, setup, savepath, magnetization_flag);
%     btpde = load_btpde(setup, savepath, magnetization_flag);
end


%% Solve BTPDE using midpoint method
if isfield(setup, "btpde_midpoint")
    disp("Computing or loading the BTPDE midpoint signals");
    savepath = create_savepath(setup, "btpde_midpoint", 'saved_simul');
    btpde_midpoint = solve_btpde_midpoint( ...
        femesh, setup, savepath, magnetization_flag ...
    );
    % btpde_midpoint = load_btpde_midpoint(setup, savepath, magnetization_flag);
end


%% Solve HADC model
if isfield(setup, "hadc")
    disp("Computing or loading the homogenized apparent diffusion coefficient");
    savepath = create_savepath(setup, "hadc");
    hadc = solve_hadc(femesh, setup, savepath);
    % hadc = load_hadc(femesh, setup, savepath);
end


%% Laplace eigendecomposition
if isfield(setup, "mf")
    % Laplace eigendecomposition
    eigenpath = create_savepath(setup, "lap_eig");
%     lap_eig = compute_laplace_eig(femesh, setup.pde, setup.mf, eigenpath);
    lap_eig = load_laplace_eig(eigenpath, setup.mf, setup.pde.mean_diffusivity);

    % Compute MF magnetization and signal
    savepath = create_savepath(setup, "mf");
    mf = solve_mf(femesh, setup, lap_eig, savepath, magnetization_flag);
%     mf = load_mf(setup, savepath, magnetization_flag);

    if isfield(setup.mf, 'hadc') && setup.mf.hadc
        mf_hadc = solve_mf_hadc(femesh, setup, lap_eig);
    end
end


%% Postprocess results
% Choose to see some of the typical plots or not
do_plots = false;

% Stop here if plotting is detoggled
if ~do_plots
    return
end

% Plot finite element mesh
if isfield(setup.geometry, "refinement")
    refinement_str = sprintf("refinement = %g", setup.geometry.refinement);
else
    refinement_str = "automatic refinement";
end
% plot_femesh(femesh, cmpts_in, cmpts_out, cmpts_ecs);
plot_femesh_everywhere(femesh, refinement_str);

% if isfield(setup, 'mf')
%     % Plot Matrix Formalism effective diffusion tensor
%     plot_diffusion_tensor(diffusion_tensor_all, setup.pde.mean_diffusivity);
% 
%     % Plot some Laplace eigenfunctions
%     if length(lap_eig) == 1
%         neig = length(lap_eig.values);
%         nshow = min(10, neig);
%         for ieig = nshow:nshow
%             diffdir = squeeze(lap_eig.moments(1, ieig, :));
%             diffdir = diffdir / norm(diffdir, 2);
%             title_str = sprintf("Laplace eigenfunction %d, l_s=%g, diffusion direction=[%.2f %.2f %.2f]",...
%                 ieig, lap_eig.length_scales(ieig), round(diffdir' * 100) / 100);
% 
%             % Split Laplace eigenfunctions into compartments
%             npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);
%             lap_eig_funcs_sep = mat2cell(lap_eig.funcs, npoint_cmpts);
%             % plot_field(femesh, lap_eig_funcs_sep, setup.pde.compartments, title_str, ieig);
%             plot_field_everywhere(femesh, lap_eig_funcs_sep, title_str, ieig);
%         end
%     else
%         % Plot the ieig-th eigenfunction for the compartment icmpt
%         icmpt = 7;
%         ncompartment = length(lap_eig);
%         icmpt = min(icmpt, ncompartment);
%         
%         ieig = 4;
%         neig = length(lap_eig(icmpt).values);
%         ieig = min(ieig, neig);
% 
%         diffdir = squeeze(lap_eig(icmpt).moments(1, ieig, :));
%         diffdir = diffdir / norm(diffdir, 2);
%         title_str = sprintf("Laplace eigenfunction %d, l_s=%g, diffusion direction=[%.2f %.2f %.2f]",...
%             icmpt, ieig, lap_eig(icmpt).length_scales(ieig), round(diffdir' * 100) / 100);
% 
%         % Split Laplace eigenfunctions into compartments
%         lap_eig_funcs_sep = cell(ncompartment, 1);
%         for jcmpt = 1:ncompartment
%             if jcmpt == icmpt
%                 lap_eig_funcs_sep{icmpt} = lap_eig(icmpt).funcs(:, ieig);
%             else
%                 lap_eig_funcs_sep{jcmpt} = lap_eig(jcmpt).funcs(:, 1) * 0;
%             end
%         end
%         plot_field_compartment(femesh, lap_eig_funcs_sep, icmpt, title_str, 1, false);
% %         plot_field(femesh, lap_eig_funcs_sep, setup.pde.compartments, title_str);
%         plot_field_everywhere(femesh, lap_eig_funcs_sep, title_str);
%     end
% end

if isfield(setup, 'mf') && isfield(setup, 'btpde')
    % Relative error between BTPDE and MF signal
    signal_allcmpts_relerr = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ...
        ./ max(abs(btpde.signal_allcmpts), [], 3);

    % Difference between BTPDE and MF signal, normalized by initial signal
    signal_allcmpts_abserr_vol = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ./ setup.pde.initial_signal;

    % Plot quantities over many directions
    if setup.ndirection > 1
        % Plot HARDI signal
        plot_hardi(setup.gradient.directions, real(btpde.signal_allcmpts) / femesh.total_volume, "BTPDE signal")
        plot_hardi(setup.gradient.directions, real(mf.signal_allcmpts) / femesh.total_volume, "MF signal")

        % Plot relative difference
        fig_title = "Rel diff between BTPDE and MF";
        plot_hardi(setup.gradient.directions, signal_allcmpts_relerr, fig_title);

        % Plot difference normalized by volume
        fig_title = "Diff between BTPDE and MF normalized by volume";
        plot_hardi(setup.gradient.directions, signal_allcmpts_abserr_vol, fig_title);
    end
end
