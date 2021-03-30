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
setup_4axons_flat;
% setup_30axons_flat;
setup_200axons;

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
ndirection = setup.gradient.ndirection;

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
if setup.gradient.flat_dirs
    flat_str = "_flat";
else
    flat_str = "";
end
dir_str = sprintf("ndir%d%s", ndirection, flat_str);
save_dir_path_spindoctor = save_dir_path_spindoctor + "/" + couple_str;

if ~isfolder(save_dir_path_spindoctor)
    mkdir(save_dir_path_spindoctor);
end


%% Solve BTPDE
if isfield(setup, "btpde")
    disp("Computing or loading the BTPDE signals");

    % Inialize BTPDE output arguments
    btpde.magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
    btpde.signal = zeros(ncompartment, namplitude, nsequence, ndirection);
    btpde.signal_allcmpts = zeros(namplitude, nsequence, ndirection);
    btpde.itertimes = zeros(namplitude, nsequence, ndirection);
    btpde.totaltime = 0;
    
    % Solve or load BTPDE, one experiment and b-value at the time
    for iseq = 1:nsequence
        for iamp = 1:namplitude
            
            % Extract experiment parameters
            seq = setup.gradient.sequences{iseq};
            bvalue = setup.gradient.bvalues(iamp, iseq);
            qvalue = setup.gradient.qvalues(iamp, iseq);

            % Create the filename string for the saving data
            if seq.delta == seq.Delta
                experi_str = sprintf("%s_dD%g", class(seq), seq.delta);
            else
                experi_str = sprintf("%s_d%g_D%g", class(seq), seq.delta, seq.Delta);
            end
            if setup.gradient.values_type == "q"
                bvalue_str = sprintf("q%g", qvalue);
            else
                bvalue_str = sprintf("b%g", bvalue);
            end
            fname = sprintf("btpde_%s_%s_%s_abstol%g_reltol%g.mat", dir_str, ...
                experi_str, bvalue_str, setup.btpde.abstol, setup.btpde.reltol);
            fname = save_dir_path_spindoctor + "/" + fname;
            
            % Choose whether to load or compute results
            if isfile(fname)
                % Load BTPDE results for one experiment and b-value
                disp("load " + fname);
                load(fname);
            else
                % Create temporary experiment structure only containing the
                % current gradient sequence and b-value
                setup_tmp = setup;
                setup_tmp.gradient.qvalues = setup.gradient.qvalues(iamp, iseq);
                setup_tmp.gradient.bvalues = setup.gradient.bvalues(iamp, iseq);
                setup_tmp.gradient.sequences = setup.gradient.sequences(iseq);

                % Solve the BTPDE for one experiment and bvalue
                btpde_tmp = solve_btpde(femesh, setup_tmp);

                % Extract results
                magnetization = btpde_tmp.magnetization;
                signal = btpde_tmp.signal;
                signal_allcmpts = btpde_tmp.signal_allcmpts;
                itertimes = btpde_tmp.itertimes;
                totaltime = btpde_tmp.totaltime;

                % Save BTPDE results
                file = save_dir_path_spindoctor + "/" + fname;
                disp("save " + fname + " -v7.3 -struct btpde_tmp");
                save(fname, "-v7.3", "-struct", "btpde_tmp");
            end

            % Store results
            btpde.magnetization(:, iamp, iseq, :) = magnetization;
            btpde.signal(:, iamp, iseq, :) = signal;
            btpde.signal_allcmpts(iamp, iseq, :) = signal_allcmpts;
            btpde.itertimes(iamp, iseq, :) = itertimes;
            btpde.totaltime = btpde.totaltime + totaltime;
        end
    end

    % BTPDE direction averaged magnetization
    btpde.magnetization_avg = average_magnetization(btpde.magnetization);

    % Clear temporary variables
    clear setup_tmp
    clear seq
    clear btpde_tmp
    clear magnetization signal signal_allcmpts itertimes totaltime
    clear fname
end


%% Solve HADC model
if isfield(setup, "hadc")
    disp("Computing or loading the homogenized apparent diffusion coefficient");

    % Initialize data structures
    hadc.adc = zeros(ncompartment, nsequence, ndirection);
    hadc.adc_allcmpts = zeros(nsequence, ndirection);
    hadc.itertimes = zeros(ncompartment, nsequence, ndirection);
    hadc.totaltime = 0;

    for iseq = 1:nsequence
        seq = setup.gradient.sequences{iseq};
        
        % File name
        fname = sprintf("hadc_%s_%s_d%g_D%g_abstol%g_reltol%g.mat", dir_str, ...
            class(seq), seq.delta, seq.Delta, setup.hadc.abstol, setup.hadc.reltol);
        fname = save_dir_path_spindoctor + "/" + fname;

        % Save or load
        if isfile(fname)
            % Load HADC results
            disp("load " + fname);
            load(fname);
        else
            % Create temporary setup structure for given index
            setup_tmp = setup;
            setup_tmp.gradient.sequences = setup.gradient.sequences(iseq);

            % Solve HADC model
            hadc_tmp = solve_hadc(femesh, setup_tmp);

            % Extract results
            adc = hadc_tmp.adc;
            adc_allcmpts = hadc_tmp.adc_allcmpts;
            itertimes = hadc_tmp.itertimes;
            totaltime = hadc_tmp.totaltime;

            % Save results
            disp("save " + fname + " -v7.3 -struct hadc_tmp");
            save(fname, "-v7.3", "-struct", "hadc_tmp");
        end

        % Store results
        hadc.adc(:, iseq, :) = adc;
        hadc.adc_allcmpts(iseq, :) = adc_allcmpts;
        hadc.itertimes(:, iseq, :) = itertimes;
        hadc.totaltime = hadc.totaltime + totaltime;
    end

    % Clear temporary variables
    clear setup_tmp
    clear seq
    clear hadc_tmp
    clear adc adc_allcmpts itertimes totaltime
    clear fname
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
        load(fname);

        % Store eigendecomposition
        lap_eig.values = values;
        lap_eig.funcs = funcs;
        lap_eig.moments = moments;
        lap_eig.massrelax = massrelax;
        lap_eig.totaltime = totaltime;
    else
        % Perform eigendecomposition
        eiglim = length2eig(setup.mf.length_scale, mean_diffusivity);
        lap_eig = compute_laplace_eig(femesh, setup.pde, eiglim, setup.mf.neig_max);

        % Save eigendecomposition
        disp("save " + fname + " -v7.3 -struct lap_eig");
        save(fname, "-v7.3", "-struct", "lap_eig");
    end

    % Clear temporary variables
    clear values funcs moments totaltime
    clear fname

    % Compute length scales of eigenvalues
    lap_eig.length_scales = eig2length(lap_eig.values, mean_diffusivity);
end


%% MF effective diffusion tensor
if isfield(setup, "mf")
    % Compute the JN value that relates the eigenmodes to their contribution
    % to the Matrix Formalism signal for a diffusion-encoding sequence
    mf_jn = compute_mf_jn(lap_eig.values, mean_diffusivity, setup);

    % Compute the Matrix Formalism effective diffusion tensor
    diffusion_tensor = compute_mf_diffusion_tensor(lap_eig, mf_jn, mean_diffusivity);
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

initial_density = setup.pde.initial_density * volumes';

% Plot finite element mesh
if isfield(setup.geometry, "refinement")
    refinement_str = sprintf("refinement = %g", setup.geometry.refinement);
else
    refinement_str = "automatic refinement";
end
% plot_femesh(femesh, cmpts_in, cmpts_out, cmpts_ecs);
plot_femesh_everywhere(femesh, refinement_str);

% Plot Matrix Formalism effective diffusion tensor
plot_diffusion_tensor(diffusion_tensor, mean_diffusivity);

% Plot some Laplace eigenfunctions
neig = length(lap_eig.values);
nshow = min(10, neig);
for ieig = nshow:nshow
    diffdir = squeeze(lap_eig.moments(1, ieig, :));
    diffdir = diffdir / norm(diffdir, 2);
    title_str = sprintf("Laplace eigenfunction %d, l_s=%g, diffusion direction=[%.2f %.2f %.2f]",...
        ieig, lap_eig.length_scales(ieig), round(diffdir' * 100) / 100);

    % Split Laplace eigenfunctions into compartments
    npoint_cmpts = cellfun(@(x) size(x, 2), femesh.points);
    lap_eig_funcs_sep = mat2cell(lap_eig.funcs, npoint_cmpts);
    % plot_field(femesh, lap_eig_funcs_sep, setup.pde.compartments, title_str, ieig);
    plot_field_everywhere(femesh, lap_eig_funcs_sep, title_str, ieig);
end

% Relative error between BTPDE and MF signal
signal_allcmpts_relerr = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ...
    ./ max(abs(btpde.signal_allcmpts), [], 3);

% Difference between BTPDE and MF signal, normalized by initial signal
signal_allcmpts_abserr_vol = abs(mf.signal_allcmpts - btpde.signal_allcmpts) ./ initial_density;

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
