%DRIVER_SPINDOCTOR Solve BTPDE, HADC, MF or analytical.
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

% Get setup
addpath setups

% setup_1axon_analytical;
% setup_1sphere_analytical;
% setup_15spheres;
% setup_2axons_deform;
% setup_5axons_myelin_relax;
% setup_neuron;
setup_4axons_flat;
% setup_30axons_flat;
% setup_30axons;
% setup_200axons;

% Choose to see some of the typical plots or not
do_plots = true;


%% Prepare experiments

% Set up the PDE model in the geometrical compartments.
setup.pde = prepare_pde(setup);

% Prepare experiments (gradient sequence, bvalues, qvalues, solvers)
setup = prepare_experiments(setup);

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

% Get sizes
ncompartment = length(setup.pde.compartments);
namplitude = length(setup.gradient.values);
nsequence = length(setup.gradient.sequences);
ndirection = setup.gradient.ndirection;


%% Perform small experiments

% Short time approximation (STA) of the ADC
[sta_adc, sta_adc_allcmpts] = compute_adc_sta(femesh, setup);

% Free diffusion signal
free = compute_free_diffusion(setup.gradient.bvalues, setup.pde.diffusivity, ...
    volumes, setup.pde.initial_density);


%% Perform analytical experiment
if isfield(setup, "analytical")
    % Solve analytical analytical model
    analytical_signal = solve_analytical(setup, volumes); % With FE volumes
    % analytical_signal = solve_analytical(setup); % With pure volumes
end


%% Solve Karger model
if isfield(setup, "karger")
    % Solve analytical analytical model
    karger = solve_karger(femesh, setup);
end


%% Perform BTPDE experiments
if isfield(setup, "btpde")
    % Solve BTPDE
    btpde = solve_btpde(femesh, setup);

    % Fit ADC from signal
    btpde_fit = fit_signal(btpde.signal, btpde.signal_allcmpts, setup.gradient.bvalues);

    % BTPDE direction averaged magnetization
    btpde.magnetization_avg = average_magnetization(btpde.magnetization);
end


%% Perform HADC experiment
if isfield(setup, "hadc")
    % Solve HADC model
    hadc = solve_hadc(femesh, setup);
end


%% Perform MF experiments
if isfield(setup, "mf")
    % Perform Laplace eigendecomposition
    eiglim = length2eig(setup.mf.length_scale, mean_diffusivity);
    lap_eig = compute_laplace_eig(femesh, setup.pde, eiglim, setup.mf.neig_max);

    % Compute length scales of eigenvalues
    lap_eig.length_scales = eig2length(lap_eig.values, mean_diffusivity);

    % Compute the JN value that relates the eigenmodes to their contribution
    % to the Matrix Formalism signal for a diffusion-encoding sequence
    mf_jn = compute_mf_jn(lap_eig.values, setup);

    % Compute the Matrix Formalism effective diffusion tensor
    diffusion_tensor = compute_mf_diffusion_tensor(femesh, lap_eig, mf_jn);

    % Compute MF magnetization
    mf = solve_mf(femesh, setup, lap_eig);

    % Fit ADC from MF signal
    mf_fit = fit_signal(mf.signal, mf.signal_allcmpts, setup.gradient.bvalues);

    % MF direction averaged magnetization
    mf.magnetization_avg = average_magnetization(mf.magnetization);
    
    % Compute MFGA signal
    mfga = compute_mfga_signal(setup, initial_signal, diffusion_tensor);
end


%% Postprocess

if ~do_plots
    return
end

if setup.geometry.cell_shape == "sphere" || setup.geometry.cell_shape == "cylinder"
    % Plot cells in canonical configuration
    plot_cells(cells, setup);
end

% Plot surface triangulation
plot_surface_triangulation(surfaces);

% Plot the finite element mesh
plot_femesh(femesh, setup.pde.compartments);
% plot_femesh_everywhere(femesh, "");

% Plot information about the geometry
plot_geometry_info(setup, volumes, surface_areas);

% Plot BTPDE magnetization in some directions
if isfield(setup, "btpde")
    for idir = 1 % :ndirection
        for iseq = 1 % :nsequence
            for iamp = 1 % :namplitude
                b = setup.gradient.bvalues(iamp, iseq);
                title_str = sprintf(...
                    "BTPDE magnetization. Sequence %d of %d, b=%.2f", ...
                    iseq, nsequence, b);
                field = btpde.magnetization(:, iamp, iseq, idir);
                % plot_field(femesh, field, setup.pde.compartments, title_str);
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
    if isfield(setup, "btpde")
        % Plot the BTPDE signal, the S0*exp(-ADC*b) curve, and the
        % free diffusion curves together.
        plot_signal(setup.gradient.bvalues, btpde.signal_allcmpts, free.signal_allcmpts, ...
            btpde_fit.S0_allcmpts, btpde_fit.adc_allcmpts, "BTPDE")

        % Plot ADC fitted from BTPDE
        plot_adc(btpde_fit.adc, btpde_fit.adc_allcmpts, "BTPDE");

        % Plot computational time
        plot_timing(btpde.itertimes, femesh, "BTPDE", "B-value");
    end

    % Plot results from HADC experiment
    if isfield(setup, "hadc")
        % Plot ADC
        plot_adc(hadc.adc, hadc.adc_allcmpts, "HADC");

        % Plot computational time
        plot_timing(hadc.itertimes, femesh, "HADC", "Compartment");
    end
else
    % Plot STA ADC
    plot_hardi(setup.gradient.directions, sta_adc_allcmpts, "STA ADC all compartments");

    % Plot BTPDE signal in all directions
    if isfield(setup, "btpde")
        % Signal before applying the gradient pulse sequence
        S0 = sum(setup.pde.initial_density .* volumes);

        % Plot normalized signals
        title_str = "BTPDE total magnetization (normalized)";
        plot_hardi(setup.gradient.directions, btpde.signal_allcmpts / S0, title_str);
    end

    % Plot HADC in all directions
    if isfield(setup, "hadc")
        % Plot normalized ADC (ADC/D0)
        title_str = sprintf("HADC all compartments");
        plot_hardi(setup.gradient.directions, hadc.adc_allcmpts / mean_diffusivity, title_str);
    end
end % One or many directions


if isfield(setup, "mf")
    % Plot Matrix Formalism effective diffusion tensor
    plot_diffusion_tensor(diffusion_tensor, mean_diffusivity);

    % Relative error between BTPDE and MF signal
    signal_allcmpts_relerr = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ...
        ./ max(abs(btpde.signal_allcmpts), [], 3);

    % Difference between BTPDE and MF signal, normalized by initial signal
    signal_allcmpts_abserr_vol = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ./ initial_signal;
    
    if ndirection == 1
        % Plot BTPDE, MF and MFGA signals
        plot_signal(setup.gradient.bvalues, btpde.signal_allcmpts, free.signal_allcmpts, btpde_fit.S0_allcmpts, btpde_fit.adc, "BTPDE signal");
        plot_signal(setup.gradient.bvalues, mf.signal_allcmpts, free.signal_allcmpts, mf_fit.S0_allcmpts, mf_fit.adc, "MF signal");
        % plot_signal(setup.gradient.bvalues, mfga.signal_allcmpts, free.signal_allcmpts, initial_signal, mfga.adc_allcmpts, "MFGA signal");
        plot_signal_btpde_mf(setup, btpde.signal_allcmpts, mf.signal_allcmpts, mfga.signal_allcmpts);

        % Display difference
        disp("Error (normalized by initial signal)");
        disp(signal_allcmpts_abserr_vol);
        disp("Relative error:");
        disp(signal_allcmpts_relerr);
    else
        % Plot HARDI signal
        plot_hardi(setup.gradient.directions, real(btpde.signal_allcmpts) / initial_signal, "BTPDE signal")
        plot_hardi(setup.gradient.directions, real(mf.signal_allcmpts) / initial_signal, "MF signal")

        % Plot relative difference
        fig_title = "Rel diff between BTPDE and MF";
        plot_hardi(setup.gradient.directions, signal_allcmpts_relerr, fig_title);

        % Plot difference normalized by volume
        fig_title = "Diff between BTPDE and MF normalized by volume";
        plot_hardi(setup.gradient.directions, signal_allcmpts_abserr_vol, fig_title);
    end
end

%%
plot_hardi(setup.gradient.directions, karger.signal_allcmpts / initial_signal, "Karger signal")
