%DRIVER_SPINDOCTOR Solve BTPDE for blood flow imaging.
clear
restoredefaultpath

% Add SpinDoctor
addpath(genpath("src"));

%% Define inputs

% Get setup
addpath(genpath("setups"));
addpath(genpath("flowterm"))

setup_1axon_flow;
%% velocity
velocity{1} = 4e-2;

%% Prepare simulation
[setup, femesh, surfaces, cells]  = prepare_simulation(setup);

% create time-independent velocity vector on each element.
velocityvec = compute_velocity(femesh, velocity, setup);

%% Perform small experiments

% Perform BTPDE midpoint experiments
if isfield(setup, "btpde_midpoint")
    % Solve BTPDE
    tautype = 1; % 0: no tau_k, 1-3: three different tau_k
    hktype = 1; % 1: directional h_k, 2:undirectional h_k
    btpde_midpoint = solve_btpde_midpoint_flow(femesh, velocityvec, tautype, hktype, setup);
end


%% Postprocess
do_plots = true;

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
plot_geometry_info(setup, femesh);


% Plot BTPDE magnetization in some directions
if isfield(setup, "btpde_midpoint")
    for idir = 1 % :setup.ndirection
        for iseq = 1 % :setup.nsequence
            for iamp = 1 % :setup.namplitude
                b = setup.gradient.bvalues(iamp, iseq);
                title_str = sprintf(...
                    "BTPDE magnetization, real part. Sequence %d of %d, b=%.2f", ...
                    iseq, setup.nsequence, b);
                field = btpde_midpoint.magnetization(:, iamp, iseq, idir);
                field{1,1} = real(field{1,1});
                % plot_field(femesh, field, setup.pde.compartments, title_str);
                plot_field_everywhere(femesh, field, title_str);
                % caxis([0 1]);
            end
        end
    end
    clear field
end

if isfield(setup, "btpde_midpoint")
    for idir = 1 % :setup.ndirection
        for iseq = 1 % :setup.nsequence
            for iamp = 1 % :setup.namplitude
                b = setup.gradient.bvalues(iamp, iseq);
                title_str = sprintf(...
                    "BTPDE magnetization, imag part. Sequence %d of %d, b=%.2f", ...
                    iseq, setup.nsequence, b);
                field = btpde_midpoint.magnetization(:, iamp, iseq, idir);
                field{1,1} = imag(field{1,1});
                % plot_field(femesh, field, setup.pde.compartments, title_str);
                plot_field_everywhere(femesh, field, title_str);
                % caxis([0 1]);
            end
        end
    end
    clear field
end
