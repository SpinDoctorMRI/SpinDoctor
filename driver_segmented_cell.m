%%DRIVER_SEGMENTED_CELL is a template for running separate simulations on
%%the dendrites and soma of a full cell.

%Add spindoctor to path
addpath(genpath('setups'));
addpath(genpath('src'));
addpath(genpath('drivers_postprocess'))
%Setup containing cell
setup_segmented_cell;

%Extract cellname and swc file
[mesh_path,cellname,~] = fileparts(setup.name);
swc_file=sprintf("swc_files/%s.swc",cellname);


%Prepare simulation
[setup, femesh, ~, ~]  = prepare_simulation(setup);
%% Segment the finite element mesh
tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
[femesh_soma,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);

ndendrites=length(femesh_dendrites);

save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);
%% Run Experiments
[sta_adc_cell, sta_adc_cell_allcmpts] = compute_adc_sta(femesh, setup);
[sta_adc_soma, sta_adc_soma_allcmpts] = compute_adc_sta(femesh, setup);
sta_adc_dendrites = cell(ndendrites,1);
for i =1:ndendrites
[sta_adc, sta_adc_allcmpts] = compute_adc_sta(femesh, setup);
sta_adc_dendrites{i} = sta_adc;
end
free_cell = compute_free_diffusion(setup.gradient.bvalues, setup.pde.diffusivity, ...
    femesh.volumes, setup.pde.initial_density);
free_soma= compute_free_diffusion(setup.gradient.bvalues, setup.pde.diffusivity, ...
    femesh_soma.volumes, setup.pde.initial_density);

free_dendrites = cell(ndendrites,1);
for i =1:ndendrites
free_dendrites{i} = compute_free_diffusion(setup.gradient.bvalues, setup.pde.diffusivity, ...
    femesh_soma.volumes, setup.pde.initial_density);
end

if isfield(setup,'btpde')
    save_path_cell = sprintf("%s/cell",save_path_root);
    % Compute BTPDE magnetization
    btpde_cell = solve_btpde(femesh, setup,save_path_cell,true);

    save_path_soma = sprintf("%s/soma",save_path_root);
    % Compute BTPDE magnetization
    btpde_soma = solve_btpde(femesh_soma, setup,save_path_soma,true);


    btpde_dendrites = cell(ndendrites,1);
    for i=1:ndendrites 
        save_path_dendrite = sprintf("%s/dendrite_%d",save_path_root,i);
        % Compute BTPDE magnetization
        btpde_dendrites{i} = solve_btpde(femesh_dendrites{i}, setup,save_path_dendrite,true);
    end
end




if isfield(setup,'mf')
    save_path_cell = sprintf("%s/cell",save_path_root);
    lap_eig = compute_laplace_eig(femesh, setup.pde, setup.mf,save_path_cell);         
    % Compute MF magnetization
    mf_cell = solve_mf_cpu(femesh, setup, lap_eig,save_path_cell,true);

    save_path_soma = sprintf("%s/soma",save_path_root);
    lap_eig = compute_laplace_eig(femesh_soma, setup.pde, setup.mf,save_path_soma);         
    % Compute MF magnetization
    mf_soma = solve_mf_cpu(femesh_soma, setup, lap_eig,save_path_soma,true);


    ndendrites=length(femesh_dendrites);
    mf_dendrites = cell(ndendrites,1);
    for i=1:ndendrites 
        save_path_dendrite = sprintf("%s/dendrite_%d",save_path_root,i);
        lap_eig = compute_laplace_eig(femesh_dendrites{i}, setup.pde, setup.mf,save_path_dendrite);         
        % Compute MF magnetization
        mf_dendrites{i} = solve_mf_cpu(femesh_dendrites{i}, setup, lap_eig,save_path_dendrite,true);
    end
end

if isfield(setup,'btpde')
    disp('BTPDE results stored in btpde_cell, btpde_soma, btpde_dendrites')
end
if isfield(setup,'mf')
    disp('MF results stored in mf_cell, mf_soma, mf_dendrites')
end
%% Post-process
addpath(genpath('drivers_postprocess'))
% Plot information about the geometry
plot_geometry_info(setup, femesh);

plot_dendrites_soma(femesh_soma,femesh_dendrites);


if isfield(setup, "btpde")
    % Fit ADC from signal
    btpde_cell_fit = fit_signal(btpde_cell.signal, btpde_cell.signal_allcmpts, setup.gradient.bvalues);
    btpde_soma_fit = fit_signal(btpde_soma.signal, btpde_soma.signal_allcmpts, setup.gradient.bvalues);
    btpde_dendrites_fit = cell(ndendrites,1);
    for i =1:ndendrites
        btpde_dendrites_fit{i} = fit_signal(btpde_dendrites{i}.signal, btpde_dendrites{i}.signal_allcmpts, setup.gradient.bvalues);
    end
    
    % Insert weighted signals
    btpde_cell.signal_weighted = btpde_cell.signal./femesh.volumes;
    btpde_soma.signal_weighted = btpde_soma.signal./femesh_soma.volumes;
    for i =1:ndendrites
        btpde_dendrites{i}.signal_weighted = btpde_dendrites{i}.signal./femesh_dendrites{i}.volumes;
    end

end


if isfield(setup, "mf")
    % Fit ADC from MF signal
    mf_cell_fit = fit_signal(mf_cell.signal, mf_cell.signal_allcmpts, setup.gradient.bvalues);
    mf_soma_fit = fit_signal(mf_soma.signal, mf_soma.signal_allcmpts, setup.gradient.bvalues);
    mf_dendrites_fit = cell(ndendrites,1);
    for i =1:ndendrites
        mf_dendrites_fit{i} = fit_signal(mf_dendrites{i}.signal, mf_dendrites{i}.signal_allcmpts, setup.gradient.bvalues);
    end

     % Insert weighted signals
    mf_cell.signal_weighted = mf_cell.signal./femesh.volumes;
    mf_soma.signal_weighted = mf_soma.signal./femesh_soma.volumes;
    for i =1:ndendrites
        mf_dendrites{i}.signal_weighted = mf_dendrites{i}.signal./femesh_dendrites{i}.volumes;
    end
end

%% Plotting results
fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
% saveas(fig,sprintf('figures/%s_tet%s.fig',cellname,setup.geometry.tetgen_options));

% Plot BTPDE magnetization in some directions
if isfield(setup, "btpde")
    for idir = 1 % :setup.ndirection
        for iseq = 1 % :setup.nsequence
            for iamp = 1 % :setup.namplitude
                b = setup.gradient.bvalues(iamp, iseq);
                % Cell
                title_str = sprintf(...
                    "BTPDE magnetization cell. Sequence %d of %d, b=%.2f", ...
                    iseq, setup.nsequence, b);
                field = btpde_cell.magnetization(:, iamp, iseq, idir);
                plot_field_everywhere(femesh, field, title_str);
                % Soma
                title_str = sprintf(...
                    "BTPDE magnetization soma. Sequence %d of %d, b=%.2f", ...
                    iseq, setup.nsequence, b);
                field = btpde_soma.magnetization(:, iamp, iseq, idir);
                plot_field_everywhere(femesh_soma, field, title_str);
                % Dendries
                for i = 1:ndendrites
                    title_str = sprintf(...
                    "BTPDE magnetization dendrite %d. Sequence %d of %d, b=%.2f", ...
                        i,iseq, setup.nsequence, b);
                    field = btpde_dendrites{i}.magnetization(:, iamp, iseq, idir);
                    plot_field_everywhere(femesh_dendrites{i}, field, title_str);
                end
            end
        end
    end
    clear field
end
if isfield(setup, "mf")
    for idir = 1 % :setup.ndirection
        for iseq = 1 % :setup.nsequence
            for iamp = 1 % :setup.namplitude
                b = setup.gradient.bvalues(iamp, iseq);
                % Cell
                title_str = sprintf(...
                    "MF magnetization cell. Sequence %d of %d, b=%.2f", ...
                    iseq, setup.nsequence, b);
                field = mf_cell.magnetization(:, iamp, iseq, idir);
                plot_field_everywhere(femesh, field, title_str);
                % Soma
                title_str = sprintf(...
                    "MF magnetization soma. Sequence %d of %d, b=%.2f", ...
                    iseq, setup.nsequence, b);
                field = mf_soma.magnetization(:, iamp, iseq, idir);
                plot_field_everywhere(femesh_soma, field, title_str);
                % Dendries
                for i = 1:ndendrites
                    title_str = sprintf(...
                    "MF magnetization dendrite %d. Sequence %d of %d, b=%.2f", ...
                        i,iseq, setup.nsequence, b);
                    field = mf_dendrites{i}.magnetization(:, iamp, iseq, idir);
                    plot_field_everywhere(femesh_dendrites{i}, field, title_str);
                end
            end
        end
    end
    clear field
end


%%
if setup.ndirection == 1
    % Plot ADC short time approximation
    plot_adc(sta_adc, sta_adc_allcmpts, "STA");

    % Plot BTPDE results
    if isfield(setup, "btpde")
        % Plot the BTPDE signal, the S0*exp(-ADC*b) curve, and the
        % free diffusion curves together.
        plot_signal(setup.gradient.bvalues, btpde_cell.signal_allcmpts, free_cell.signal_allcmpts, ...
            btpde_cell_fit.S0_allcmpts, btpde_cell_fit.adc_allcmpts, "BTPDE cell")

        % Plot ADC fitted from BTPDE
        plot_adc(btpde_cell_fit.adc, btpde_cell_fit.adc_allcmpts, "BTPDE cell");

        % Plot computational time
        plot_timing(btpde_cell.itertimes, femesh, "BTPDE cell", "B-value");

        % Plot the BTPDE signal, the S0*exp(-ADC*b) curve, and the
        % free diffusion curves together.
        plot_signal(setup.gradient.bvalues, btpde_soma.signal_allcmpts, free_soma.signal_allcmpts, ...
            btpde_soma_fit.S0_allcmpts, btpde_soma_fit.adc_allcmpts, "BTPDE soma")

        % Plot ADC fitted from BTPDE
        plot_adc(btpde_soma_fit.adc, btpde_soma_fit.adc_allcmpts, "BTPDE soma");

        % Plot computational time
        plot_timing(btpde_soma.itertimes, femesh_soma, "BTPDE soma", "B-value");

        for i =1:ndendrites
        % Plot the BTPDE signal, the S0*exp(-ADC*b) curve, and the
        % free diffusion curves together.
        plot_signal(setup.gradient.bvalues, btpde_dendrites{i}.signal_allcmpts, free_dendrites{i}.signal_allcmpts, ...
            btpde_dendrites_fit{i}.S0_allcmpts, btpde_dendrites_fit{i}.adc_allcmpts, "BTPDE dendrite "+string(i))

        % Plot ADC fitted from BTPDE
        plot_adc(btpde_dendrites_fit{i}.adc, btpde_dendrites_fit{i}.adc_allcmpts,  "BTPDE dendrite "+string(i));

        % Plot computational time
        plot_timing(btpde_dendrites{i}.itertimes, femesh_dendrites{i},  "BTPDE dendrite "+string(i), "B-value");
        end
    end
else
    % Plot STA ADC
    

    % Plot BTPDE signal in all directions
    if isfield(setup, "btpde")
        % Plot normalized signals
        plot_hardi(setup.gradient.directions, sta_adc_cell_allcmpts, "STA ADC cell");
        title_str = "BTPDE cell total magnetization (normalized)";
        plot_hardi(setup.gradient.directions, btpde_cell.signal_weighted, title_str);
        plot_hardi(setup.gradient.directions, sta_adc_soma_allcmpts, "STA ADC soma");
        title_str = "BTPDE soma total magnetization (normalized)";
        plot_hardi(setup.gradient.directions, btpde_soma.signal_weighted, title_str);
        for i =1:ndendrites
        plot_hardi(setup.gradient.directions, sta_adc_dendrites{i}, "STA ADC dendrite "+string(i));
        title_str = sprintf("BTPDE dendrite %d total magnetization (normalized)",i);
        plot_hardi(setup.gradient.directions, btpde_dendrites{i}.signal_weighted , title_str);
        end
    end

end % One or many directions

