%DRIVER_SPINDOCTOR Solve BTPDE, HADC, MF or multilayer.
%   Compare different ADC. Plot results in many directions.
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

% Define the directory where the input file is kept
addpath input_files

% Load input parameters using input file script. The following
% structures will be added to the workspace:
%    params_cells
%    params_domain
%    experiment
% Choose input script:

% define_params_1axon_multilayer;
% define_params_1sphere_multilayer;
% define_params_15spheres;
% define_params_2axons_deform;
% define_params_5axons_myelin_relax;
% define_params_neuron;
define_params_4axons_flat;
% define_params_30axons_flat;

% Choose to see some of the typical plots or not
do_plots = true;


%% Prepare experiments

% Create geometrical configuration
cells = create_cells(params_cells);

% Set up the PDE model in the geometrical compartments.
params_domain = prepare_pde(params_cells, params_domain);

% Prepare experiments (gradient sequence, bvalues, qvalues)
experiment = prepare_experiments(experiment);

% Create or load finite element mesh
[femesh, surfaces] = create_femesh(cells, params_cells, params_domain);

% Get volume and surface area quantities from mesh
[volumes, surface_areas] = get_vol_sa(femesh);

% Compute volume weighted mean of diffusivities over compartments
% Take trace of each diffusion tensor, divide by 3
mean_diffusivity = trace(sum(params_domain.diffusivity .* shiftdim(volumes, -1), 3)) ...
    / (3 * sum(volumes));

% Initial total signal
initial_signal = params_domain.initial_density * volumes';


%% Perform experiments

% Short time approximation (STA) of the ADC
[sta_adc, sta_adc_allcmpts] = compute_adc_sta(femesh, params_domain, experiment);

% Free diffusion signal
free = compute_free_diffusion(experiment.bvalues, params_domain.diffusivity, ...
    volumes, params_domain.initial_density);

% Perform BTPDE experiments
if isfield(experiment, "btpde")
    % Solve BTPDE
    btpde = solve_btpde(femesh, params_domain, experiment);

    % Fit ADC from signal
    btpde_fit = fit_signal(btpde.signal, btpde.signal_allcmpts, experiment.bvalues);

    % BTPDE direction averaged magnetization
    btpde.magnetization_avg = average_magnetization(btpde.magnetization);
end

% Perform HADC experiment
if isfield(experiment, "hadc")
    % Solve HADC model
    hadc = solve_hadc(femesh, params_domain, experiment);
end

% Perform MF experiments
if isfield(experiment, "mf")
    % Perform Laplace eigendecomposition
    eiglim = length2eig(experiment.mf.length_scale, mean_diffusivity);
    lap_eig = compute_laplace_eig(femesh, params_domain, eiglim, experiment.mf.neig_max);

    % Compute length scales of eigenvalues
    lap_eig.length_scales = eig2length(lap_eig.values, mean_diffusivity);

    % Compute the JN value that relates the eigenmodes to their contribution
    % to the Matrix Formalism signal for a diffusion-encoding sequence
    mf_jn = compute_mf_jn(lap_eig.values, mean_diffusivity, experiment);

    % Compute the Matrix Formalism effective diffusion tensor
    diffusion_tensor = compute_mf_diffusion_tensor(lap_eig, mf_jn, mean_diffusivity);

    % Compute MF magnetization
    mf = solve_mf(femesh, params_domain, experiment, lap_eig);

    % Fit ADC from MF signal
    mf_fit = fit_signal(mf.signal, mf.signal_allcmpts, experiment.bvalues);

    % MF direction averaged magnetization
    mf.magnetization_avg = average_magnetization(mf.magnetization);
    
    % Compute MFGA signal
    mfga = compute_mfga_signal(experiment, initial_signal, diffusion_tensor);
end


%% Postprocess

if ~do_plots
    return
end

if params_cells.shape == "sphere" || params_cells.shape == "cylinder"
    % Plot cells in canonical configuration
    plot_cells(cells, params_cells);
end

% Plot surface triangulation
plot_surface_triangulation(surfaces);

% Plot the finite element mesh
plot_femesh(femesh, params_domain.compartments);
% plot_femesh_everywhere(femesh, "");

% Plot information about the geometry
plot_geometry_info(params_domain, volumes, surface_areas);

% Get sizes
ncompartment = length(params_domain.compartments);
namplitude = length(experiment.values);
nsequence = length(experiment.sequences);
ndirection = experiment.ndirection;

% Plot BTPDE magnetization in some directions
if isfield(experiment, "btpde")
    for idir = 1 % :ndirection
        for iseq = 1 % :nsequence
            for iamp = 1 % :namplitude
                b = experiment.bvalues(iamp, iseq);
                title_str = sprintf(...
                    "BTPDE magnetization. Sequence %d of %d, b=%.2f", ...
                    iseq, nsequence, b);
                field = btpde.magnetization(:, iamp, iseq, idir);
                % plot_field(femesh, field, params_domain.compartments, title_str);
                plot_field_everywhere(femesh, field, title_str);
                % caxis([0 1]);
            end
        end
    end
    clear field
end

% Plot results in one or many directions
if ndirection == 1
    % Plot ADC short time approximation
    plot_adc(sta_adc, sta_adc_allcmpts, "STA");

    % Plot BTPDE results
    if isfield(experiment, "btpde")
        % Plot the BTPDE signal, the S0*exp(-ADC*b) curve, and the
        % free diffusion curves together.
        plot_signal(experiment.bvalues, btpde.signal_allcmpts, free.signal_allcmpts, ...
            btpde_fit.S0_allcmpts, btpde_fit.adc_allcmpts, "BTPDE")

        % Plot ADC fitted from BTPDE
        plot_adc(btpde_fit.adc, btpde_fit.adc_allcmpts, "BTPDE");

        % Plot computational time
        plot_timing(btpde.itertimes, femesh, "BTPDE", "B-value");
    end

    % Plot results from HADC experiment
    if isfield(experiment, "hadc")
        % Plot ADC
        plot_adc(hadc.adc, hadc.adc_allcmpts, "HADC");

        % Plot computational time
        plot_timing(hadc.itertimes, femesh, "HADC", "Compartment");
    end
else
    % Plot STA ADC
    plot_hardi(experiment.directions, sta_adc_allcmpts, "STA ADC all compartments");

    % Plot BTPDE signal in all directions
    if isfield(experiment, "btpde")
        % Signal before applying the gradient pulse sequence
        S0 = sum(params_domain.initial_density .* volumes);

        % Plot normalized signals
        title_str = "BTPDE total magnetization (normalized)";
        plot_hardi(experiment.directions, btpde.signal_allcmpts / S0, title_str);
    end

    % Plot HADC in all directions
    if isfield(experiment, "hadc")
        % Plot normalized ADC (ADC/D0)
        title_str = sprintf("HADC all compartments");
        plot_hardi(experiment.directions, hadc.adc_allcmpts / mean_diffusivity, title_str);
    end
end % One or many directions


if isfield(experiment, "mf")
    % Plot Matrix Formalism effective diffusion tensor
    plot_diffusion_tensor(diffusion_tensor, mean_diffusivity);

    % Relative error between BTPDE and MF signal
    signal_allcmpts_relerr = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ...
        ./ max(abs(btpde.signal_allcmpts), [], 3);

    % Difference between BTPDE and MF signal, normalized by initial signal
    signal_allcmpts_abserr_vol = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ./ initial_signal;
    
    if ndirection == 1
        % Plot BTPDE, MF and MFGA signals
        plot_signal(experiment.bvalues, btpde.signal_allcmpts, free.signal_allcmpts, btpde_fit.S0_allcmpts, btpde_fit.adc, "BTPDE signal");
        plot_signal(experiment.bvalues, mf.signal_allcmpts, free.signal_allcmpts, mf_fit.S0_allcmpts, mf_fit.adc, "MF signal");
        % plot_signal(experiment.bvalues, mfga.signal_allcmpts, free.signal_allcmpts, initial_signal, mfga.adc_allcmpts, "MFGA signal");
        plot_signal_btpde_mf(experiment, btpde.signal_allcmpts, mf.signal_allcmpts, mfga.signal_allcmpts);

        % Display difference
        disp("Error (normalized by initial signal)");
        disp(signal_allcmpts_abserr_vol);
        disp("Relative error:");
        disp(signal_allcmpts_relerr);
    else
        % Plot HARDI signal
        plot_hardi(experiment.directions, real(btpde.signal_allcmpts) / initial_signal, "BTPDE signal")
        plot_hardi(experiment.directions, real(mf.signal_allcmpts) / initial_signal, "MF signal")

        % Plot relative difference
        fig_title = "Rel diff between BTPDE and MF";
        plot_hardi(experiment.directions, signal_allcmpts_relerr, fig_title);

        % Plot difference normalized by volume
        fig_title = "Diff between BTPDE and MF normalized by volume";
        plot_hardi(experiment.directions, signal_allcmpts_abserr_vol, fig_title);
    end
end
