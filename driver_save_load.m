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
define_params_2axons_deform;
% define_params_5axons_myelin_relax;
% define_params_neuron;
% define_params_4axons_flat;
% define_params_30axons_flat;

% Choose to see some of the typical plots or not
do_plots = true;


%% Prepare experiments

% Create the geometrical configuration
cells = create_cells(params_cells);

% Set up the PDE model in the geometrical compartments
params_domain = prepare_pde(params_cells, params_domain);

% Prepare experiments (gradient sequence, bvalues, qvalues)
experiment = prepare_experiments(experiment);

% Get sizes
ncompartment = length(params_domain.compartments);
namplitude = length(experiment.values);
nsequence = length(experiment.sequences);
ndirection = experiment.ndirection;

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


% Create folder for saving results
tmp = split(params_cells.filename, "/");
tmp = tmp(end);
if endsWith(tmp, ".stl")
    tmp = split(tmp, ".stl");
    tmp = tmp(1);
end
refinement_str = "";
if isfield(params_cells, "refinement")
    refinement_str = sprintf("_refinement%g", params_cells.refinement);
end
save_dir_path_spindoctor = "saved_simul/" + tmp + refinement_str;
couple_str = sprintf("kappa%g_%g", params_domain.permeability_in_out, ...
    params_domain.permeability_out_ecs);
if experiment.flat_dirs
    flat_str = "_flat";
else
    flat_str = "";
end
save_dir_path_spindoctor = sprintf("%s/%s_ndir%d%s", save_dir_path_spindoctor, ...
    couple_str, ndirection, flat_str);

if ~isfolder(save_dir_path_spindoctor)
    mkdir(save_dir_path_spindoctor);
end


%% Solve BTPDE
if isfield(experiment, "btpde")
    disp("Computing or loading the BTPDE signals");

    % Inialize BTPDE output arguments
    btpde.magnetization = cell(ncompartment, namplitude, nsequence, ndirection);
    btpde.signal = zeros(ncompartment, namplitude, nsequence, ndirection);
    btpde.signal_allcmpts = zeros(namplitude, nsequence, ndirection);
    btpde.itertimes = zeros(namplitude, nsequence, ndirection);
    btpde.totaltime = 0;

    % Solve or load BTPDE, one experiment and b-value at the time
    for iamp = 1:namplitude
        for iseq = 1:nsequence
            % Extract experiment parameters
            seq = experiment.sequences{iseq};
            bvalue = experiment.bvalues(iamp, iseq);
            qvalue = experiment.qvalues(iamp, iseq);

            % Create the filename string for the saving data
            if seq.delta == seq.Delta
                experi_str = sprintf("%s_dD%g", class(seq), seq.delta);
            else
                experi_str = sprintf("%s_d%g_D%g", class(seq), seq.delta, seq.Delta);
            end
            if experiment.values_type == "q"
                bvalue_str = sprintf("q%g", qvalue);
            else
                bvalue_str = sprintf("b%g", bvalue);
            end
            fname = sprintf("btpde_%s_%s_abstol%g_reltol%g.mat", experi_str, ...
                bvalue_str, experiment.btpde.abstol, experiment.btpde.reltol);
            fname = save_dir_path_spindoctor + "/" + fname;
            
            % Choose whether to load or compute results
            if isfile(fname)
                % Load BTPDE results for one experiment and b-value
                disp("load " + fname);
                load(fname);
            else
                % Create temporary experiment structure only containing the
                % current gradient sequence and b-value
                experi_tmp = experiment;
                experi_tmp.qvalues = experiment.qvalues(iamp, iseq);
                experi_tmp.bvalues = experiment.bvalues(iamp, iseq);
                experi_tmp.sequences = experiment.sequences(iseq);

                % Solve the BTPDE for one experiment and bvalue
                btpde_tmp = solve_btpde(femesh, params_domain, experi_tmp);

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
    clear experi_tmp
    clear seq
    clear btpde_tmp
    clear magnetization signal signal_allcmpts itertimes totaltime
    clear fname
end


%% Solve HADC model
if isfield(experiment, "hadc")
    disp("Computing or loading the homogenized apparent diffusion coefficient");

    % Initialize data structures
    hadc.adc = zeros(ncompartment, nsequence, ndirection);
    hadc.adc_allcmpts = zeros(nsequence, ndirection);
    hadc.itertimes = zeros(ncompartment, nsequence, ndirection);
    hadc.totaltime = 0;

    for iseq = 1:nsequence
        seq = experiment.sequences{iseq};
        
        % File name
        fname = sprintf("hadc_%s_d%g_D%g_abstol%g_reltol%g.mat",...
            class(seq), seq.delta, seq.Delta, experiment.hadc.abstol, experiment.hadc.reltol);
        fname = save_dir_path_spindoctor + "/" + fname;

        % Save or load
        if isfile(fname)
            % Load HADC results
            disp("load " + fname);
            load(fname);
        else
            % Create temporary experiment structure for given experiment index
            experi_tmp = experiment;
            experi_tmp.sequences = experiment.sequences(iseq);

            % Solve HADC model
            hadc_tmp = solve_hadc(femesh, params_domain, experi_tmp);

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
        hadc.totaltime(:, iseq, :) = hadc.totaltime(:, iseq, :) + totaltime;
    end

    % Clear temporary variables
    clear experi_tmp
    clear seq
    clear hadc_tmp
    clear adc adc_allcmpts itertimes totaltime
    clear fname
end


%% Laplace eigendecomposition
if isfield(experiment, "mf")
    disp("Computing or loading the Laplace eigenfunctions");

    % Filename
    fname = sprintf("lap_eig_lengthscale%g.mat", experiment.mf.length_scale);
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
        lap_eig.totaltime = totaltime;
    else
        % Perform eigendecomposition
        eiglim = length2eig(experiment.mf.length_scale, mean_diffusivity);
        lap_eig = compute_laplace_eig(femesh, params_domain, eiglim, experiment.mf.neig_max);

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
if isfield(experiment, "mf")
    % Compute the JN value that relates the eigenmodes to their contribution
    % to the Matrix Formalism signal for a diffusion-encoding sequence
    mf_jn = compute_mf_jn(lap_eig.values, mean_diffusivity, experiment);

    % Compute the Matrix Formalism effective diffusion tensor
    diffusion_tensor = compute_mf_diffusion_tensor(lap_eig, mf_jn, mean_diffusivity);
end


%% Compute MF magnetization
if isfield(experiment, "mf")
    % Compute MF magnetization and signal
    mf = solve_mf(femesh, params_domain, experiment, lap_eig);

    % MF direction averaged magnetization
    mf.magnetization_avg = average_magnetization(mf.magnetization);
end


%% Postprocess results

% Stop here if plotting is detoggled
if ~do_plots
    return
end

initial_density = params_domain.initial_density * volumes';

% Plot finite element mesh
% plot_femesh(femesh, cmpts_in, cmpts_out, cmpts_ecs);
if isfield(params_cells, "refinement")
    refinement_str = sprintf("refinement = %g", params_domain.refinement);
else
    refinement_str = "automatic refinement";
end
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
    % plot_field(femesh, lap_eig_funcs_sep, cmpts_in, cmpts_out, cmpts_ecs, title_str, ieig);
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
    plot_hardi(experiment.directions, real(btpde.signal_allcmpts) / sum(volumes), "BTPDE signal")
    plot_hardi(experiment.directions, real(mf.signal_allcmpts) / sum(volumes), "MF signal")

    % Plot relative difference
    fig_title = "Rel diff between BTPDE and MF";
    plot_hardi(experiment.directions, signal_allcmpts_relerr, fig_title);

    % Plot difference normalized by volume
    fig_title = "Diff between BTPDE and MF normalized by volume";
    plot_hardi(experiment.directions, signal_allcmpts_abserr_vol, fig_title);
end
