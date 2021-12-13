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

setup_4axons_flat_vargdir;


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
ndirection = size(setup.gradient.directions, 2);


%% Perform small experiments

% % Short time approximation (STA) of the ADC
% [sta_adc, sta_adc_allcmpts] = compute_adc_sta(femesh, setup);

% Free diffusion signal
free = compute_free_diffusion(setup.gradient.bvalues, setup.pde.diffusivity, ...
    volumes, setup.pde.initial_density);


%% Perform analytical experiment
if isfield(setup, "analytical")
    % Solve analytical analytical model
    analytical_signal = solve_analytical(setup, volumes); % With FE volumes
    % analytical_signal = solve_analytical(setup); % With pure volumes
end





%% Perform BTPDE experiments
if isfield(setup, "btpde")
    % Solve BTPDE
    
    
    %%%%%%%%%% For CoreyBaron %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% loads the gradient steps
    load gradientWaveforms.mat
    
    %%% rotate the geometry to check the signal stays the same
    rotationangle = 0.35*pi;
    for icmpt = 1:ncompartment
        points_old = femesh.points{icmpt};
        points = points_old;
        points(1,:) = points_old(1,:)*cos(rotationangle) + points_old(2,:)*(-sin(rotationangle));
        points(2,:) = points_old(1,:)*sin(rotationangle) + points_old(2,:)*(cos(rotationangle));
        femesh.points{icmpt} = points;
    end
    
    %%% this is the list of times, starting at 0, every 20 milliseconds
    timelist_vargdir = 0:20:20*length(grads);
    
    %%% convert gradient to q = gamma*g
    vargdir = grads*2.67513e-7*1000;
    
    %%% add the time steps and the gradients to setup
    setup.vargdir = vargdir;
    setup.timelist_vargdir = timelist_vargdir;
    
    %%% modified solve_btpde function to allow var gdir and called it
    %%% solve_btpde_vargdir
    btpde = solve_btpde_vargdir(femesh,setup);
    
    %%%%%%%%%%%%% End of change for CoreyBaron %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Fit ADC from signal
    btpde_fit = fit_signal(btpde.signal, btpde.signal_allcmpts, setup.gradient.bvalues);
    
    % BTPDE direction averaged magnetization
    btpde.magnetization_avg = average_magnetization(btpde.magnetization);
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
%plot_femesh(femesh, setup.pde.compartments);
plot_femesh_everywhere(femesh, "");

% Plot information about the geometry
plot_geometry_info(setup, femesh);

% Plot BTPDE magnetization in some directions
if isfield(setup, "btpde")
    for idir = 1 % :ndirection
        for iseq = 1:nsequence
            for iamp = 1:namplitude
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
    %     plot_adc(sta_adc, sta_adc_allcmpts, "STA");
    
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


