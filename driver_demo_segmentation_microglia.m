clear;
clc;
tic
addpath(genpath('setups'));
addpath(genpath('src'));
addpath(genpath('drivers_postprocess'));


setup_file = "setup_segmentation_microglia";


% Launch setup here
run(sprintf("%s.m",setup_file));
fprintf("Running %s.m\n",setup_file)
%% Obtain finite element meshes

[~,cellname,~] = fileparts(setup.name);


[setup, femesh_cell, ~, ~,femesh_soma,femesh_neurites]  = prepare_simulation(setup);

%% Run simulations


save_magnetization = true; include_cell = true;

[btpde_cell,btpde_soma,btpde_neurites] = run_btpde_cell(femesh_cell, setup,save_magnetization,femesh_soma,femesh_neurites,include_cell);
[mf_cell,mf_soma,mf_neurites,lap_eig_cell,lap_eig_soma,lap_eig_neurites] = run_mf_cell(femesh_cell, setup,save_magnetization,femesh_soma,femesh_neurites,include_cell);

%% Load simulations
% 
load_magnetization = true; include_cell = true; include_soma= true;
nneurites= length(femesh_neurites);

[btpde_cell,btpde_soma,btpde_neurites] =load_btpde_cell(setup,load_magnetization,nneurites,include_soma,include_cell);
[mf_cell,mf_soma,mf_neurites] = load_btpde_cell(setup,load_magnetization,nneurites,include_soma,include_cell);


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


%%Version which is single line for batch processing
mesh="mesh_files/microglia/Ramified/826_5_3.ply";
setup_file="setup_camino"; %% Change setup file to change sequences and parameters of mf/btpde, permeability, diffusivity etc.
tetgen_options="-pq1.2a1.0O9VCn";
swc_file="swc_files/826_5_3.swc";
soma_file="mesh_files/microglia/Soma/826_5_3.ply";
[results,femesh_cell,femesh_soma,femesh_neurites]= run_simulations_microglia(mesh,setup_file,tetgen_options,swc_file,soma_file);

%%Load results
[results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_microglia(mesh,setup_file,tetgen_options,swc_file,soma_file);