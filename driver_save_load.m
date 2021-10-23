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
save_magnetization = false;

% Choose to see some of the typical plots or not
do_plots = true;


%% Prepare experiments

% Set up the PDE model in the geometrical compartments
setup.pde = prepare_pde(setup);

% Prepare experiments (gradient sequence, bvalues, qvalues)
setup = prepare_experiments(setup);

% Get sizes
ncompartment = length(setup.pde.compartments);
namplitude = length(setup.gradient.values);
nsequence = length(setup.gradient.sequences);
ndirection = size(setup.gradient.directions, 2);

% Create or load finite element mesh
[femesh, surfaces, cells] = create_geometry(setup);

% Get volume and surface area quantities from mesh
[volumes, surface_areas] = get_vol_sa(femesh);

% Compute volume weighted mean of diffusivities over compartments
% Take trace of each diffusion tensor, divide by 3
mean_diffusivity = trace(sum(setup.pde.diffusivity .* shiftdim(volumes, -1), 3)) ...
    / (3 * sum(volumes));

% Initial total signal
initial_signal = setup.pde.initial_density * volumes';

% Create folder for saving results
tmp = split(setup.name, "/");
tmp = tmp(end);
if endsWith(tmp, ".stl")
    tmp = split(tmp, ".stl");
    tmp = tmp(1);
end
refinement_str = "";
if isfield(setup.geometry, "refinement")
    refinement_str = sprintf("_refinement%g", setup.geometry.refinement);
end
save_dir_path_spindoctor = "saved_simul/" + tmp + refinement_str;
couple_str = sprintf("kappa%g_%g", setup.pde.permeability_in_out, ...
    setup.pde.permeability_out_ecs);
dir_str = sprintf("ndir%d", ndirection);
save_dir_path_spindoctor = save_dir_path_spindoctor + "/" + couple_str;

if ~isfolder(save_dir_path_spindoctor)
    mkdir(save_dir_path_spindoctor);
end


%% Solve BTPDE
if isfield(setup, "btpde")
    disp("Computing or loading the BTPDE signals");
    btpde = solve_btpde(femesh, setup, save_dir_path_spindoctor, save_magnetization);
    btpde.magnetization_avg = average_magnetization(btpde.magnetization);
end


%% Solve BTPDE using midpoint method
if isfield(setup, "btpde_midpoint")
    disp("Computing or loading the BTPDE midpoint signals");
    btpde_midpoint = solve_btpde_midpoint( ...
        femesh, setup, save_dir_path_spindoctor, save_magnetization ...
    );
    btpde_midpoint.magnetization_avg ...
        = average_magnetization(btpde_midpoint.magnetization);
end


%% Solve HADC model
if isfield(setup, "hadc")
    disp("Computing or loading the homogenized apparent diffusion coefficient");
    hadc = solve_hadc(femesh, setup, save_dir_path_spindoctor);
end


%% Laplace eigendecomposition
if isfield(setup, "mf")
    disp("Computing or loading the Laplace eigenfunctions");

    % Filename
    fname = sprintf("lap_eig_lengthscale%g.mat", setup.mf.length_scale);
    fname = save_dir_path_spindoctor + "/" + fname;
    
    % Save or load
    if isfile(fname)
        % Load eigendecomposition
        disp("load " + fname);
        lap_eig = load(fname);

    else
        % Perform eigendecomposition
        eiglim = length2eig(setup.mf.length_scale, mean_diffusivity);
        lap_eig = compute_laplace_eig(femesh, setup.pde, setup.mf, eiglim);

        % Save eigendecomposition
        disp("save " + fname + " -v7.3 -struct lap_eig");
        save(fname, "-v7.3", "-struct", "lap_eig");
    end

    % Clear temporary variables
    clear values funcs moments totaltime
    clear fname

    % Compute length scales of eigenvalues
    lap_eig = add_eig_length(lap_eig, mean_diffusivity);
end


%% MF effective diffusion tensor
if isfield(setup, "mf")
    % Compute the JN value that relates the eigenmodes to their contribution
    % to the Matrix Formalism signal for a diffusion-encoding sequence
    mf_jn = compute_mf_jn(lap_eig, setup);

    % Compute the Matrix Formalism effective diffusion tensor
    [diffusion_tensor, diffusion_tensor_all] = compute_mf_diffusion_tensor(femesh, lap_eig, mf_jn);
end


%% Compute MF magnetization
if isfield(setup, "mf")
    % Compute MF magnetization and signal
    mf = solve_mf(femesh, setup, lap_eig);

    % MF direction averaged magnetization
    mf.magnetization_avg = average_magnetization(mf.magnetization);
end


%% Postprocess results

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

% Plot Matrix Formalism effective diffusion tensor
plot_diffusion_tensor(diffusion_tensor_all, mean_diffusivity);

% Plot some Laplace eigenfunctions
% TO DO: add support for uncorrelated compartments
% neig = length(lap_eig.values);
% nshow = min(10, neig);
% for ieig = nshow:nshow
%     diffdir = squeeze(lap_eig.moments(1, ieig, :));
%     diffdir = diffdir / norm(diffdir, 2);
%     title_str = sprintf("Laplace eigenfunction %d, l_s=%g, diffusion direction=[%.2f %.2f %.2f]",...
%         ieig, lap_eig.length_scales(ieig), round(diffdir' * 100) / 100);

%     % Split Laplace eigenfunctions into compartments
%     npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);
%     lap_eig_funcs_sep = mat2cell(lap_eig.funcs, npoint_cmpts);
%     % plot_field(femesh, lap_eig_funcs_sep, setup.pde.compartments, title_str, ieig);
%     plot_field_everywhere(femesh, lap_eig_funcs_sep, title_str, ieig);
% end

% Relative error between BTPDE and MF signal
signal_allcmpts_relerr = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ...
    ./ max(abs(btpde.signal_allcmpts), [], 3);

% Difference between BTPDE and MF signal, normalized by initial signal
signal_allcmpts_abserr_vol = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ./ initial_signal;

% Plot quantities over many directions
if ndirection > 1
    % Plot HARDI signal
    plot_hardi(setup.gradient.directions, real(btpde.signal_allcmpts) / sum(volumes), "BTPDE signal")
    plot_hardi(setup.gradient.directions, real(mf.signal_allcmpts) / sum(volumes), "MF signal")

    % Plot relative difference
    fig_title = "Rel diff between BTPDE and MF";
    plot_hardi(setup.gradient.directions, signal_allcmpts_relerr, fig_title);

    % Plot difference normalized by volume
    fig_title = "Diff between BTPDE and MF normalized by volume";
    plot_hardi(setup.gradient.directions, signal_allcmpts_abserr_vol, fig_title);
end
