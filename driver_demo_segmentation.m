clear;
clc;
tic
addpath(genpath('setups'));
addpath(genpath('src'));
addpath(genpath('drivers_postprocess'));


setup_file = "setup_segmentation";


% Launch setup here
run(sprintf("%s.m",setup_file));
fprintf("Running %s.m\n",setup_file)
%% Obtain finite element meshes

[~,cellname,~] = fileparts(setup.name);

[setup, femesh_cell, ~, ~,femesh_soma,femesh_neurites]  = prepare_simulation(setup);

%% Run simulations
savepath_root=sprintf("saved_simul\\%s_tet%s",cellname,setup.geometry.tetgen_options);

save_magnetization = true; include_cell = true;

[btpde_cell,btpde_soma,btpde_neurites] = run_btpde_cell(femesh_cell, setup, savepath_root,save_magnetization,femesh_soma,femesh_neurites,include_cell);
[mf_cell,mf_soma,mf_neurites,lap_eig_cell,lap_eig_soma,lap_eig_neurites] = run_mf_cell(femesh_cell, setup, savepath_root,save_magnetization,femesh_soma,femesh_neurites,include_cell);

%% Visualize results



plot_neurites_soma(femesh_soma,femesh_neurites);


ifield = 1;

plot_field_everywhere(femesh_cell, btpde_cell.magnetization, sprintf('Cell BTPDE Magentization, b = %.1f',setup.gradient.bvalues(ifield)), ifield);
plot_field_everywhere(femesh_cell, mf_cell.magnetization, sprintf('Cell MF Magentization, b = %.1f',setup.gradient.bvalues(ifield)), ifield);


plot_field_everywhere(femesh_soma, btpde_soma.magnetization, sprintf('Soma BTPDE Magentization, b = %.1f',setup.gradient.bvalues(ifield)), ifield);
plot_field_everywhere(femesh_soma, btpde_soma.magnetization, sprintf('Soma MF Magentization, b = %.1f',setup.gradient.bvalues(ifield)), ifield);

nneruites=length(femesh_neurites);
for ib = 1:nneruites
plot_field_everywhere(femesh_neurites{ib}, btpde_neurites{ib}.magnetization, sprintf('Neurite %d BTPDE Magentization, b = %.1f',ib,setup.gradient.bvalues(ifield)), ifield);
plot_field_everywhere(femesh_neurites{ib}, mf_neurites{ib}.magnetization, sprintf('Neurite %d MF Magentization, b = %.1f',ib,setup.gradient.bvalues(ifield)), ifield);
end