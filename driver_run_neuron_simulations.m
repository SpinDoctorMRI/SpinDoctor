%% Add SpinDoctor to Path
restoredefaultpath;
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'))

setup_file='setup_pgse_neuron'; tetgen_options = "-pq1.2a0.5O9VCn";
mesh='mesh_files/selected/1-2-2.CNG.ply';
[~,cellname,~] = fileparts(mesh);

mesh_um = sprintf("mesh_files/ultraliser_modified/%s_um.ply",cellname);
%% Run simulations 

% Alpha_Mesh_Swc mesh
[results,femesh_cell,~,~]=run_simulations_neuron(mesh,setup_file,tetgen_options);

% Modified Ultraliser mesh
[results_um,femesh_cell_um,~,~]= run_simulations_neuron(mesh_um,setup_file,tetgen_options);

%% Load simulations
% [results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
% [results_um,femesh_cell_um,~,~]= load_simulations_neuron(mesh_um,setup_file,tetgen_options);
%% Visualise and compare signals

% Extract bvalues and signals
bvals = results.setup.gradient.bvalues;
signal = real(results.mf_cell.signal./femesh_cell.total_volume);
signal_um = real(results_um.mf_cell.signal./femesh_cell_um.total_volume);


volume=femesh_cell.total_volume;
volume_um = femesh_cell_um.total_volume;

filename=sprintf('neuron_meshing_paper/neuron_output/%s.mat',cellname);
fprintf('Saving signals and volume to %s\n',filename);
save(filename,'signal','signal_um','bvals','volume','volume_um');
namplitude = length(bvals);

rel_diff = abs((signal- signal_um)./signal);
abs_diff= abs(signal- signal_um);

figure;

figure;
plot(bvals,signal);
hold on;
plot(bvals,signal_um);
title("Comparing volume weighted signals");
ylabel("Volume weighted signals");
legend(["Alpha_Mesh_Swc","Modified Ultraliser"],'Interpreter','none');
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');


figure;
plot(bvals,rel_diff);
title("Absolute relative difference");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

figure;
plot(bvals,abs_diff);
title("Absolute difference");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

fprintf("Alpha_Mesh_Swc volume = %f,\n Modified Ultraliser volume = %f.\n",volume,volume_um);
fprintf("Here we only plot and consider the volume-weighted signals:\n")
fprintf("Maximum relative signal difference = %f,\n Maximum absolute signal difference =%f.\n",max(rel_diff),max(abs_diff));


%% Run segmented cell simulations
tetgen_options = "-pq1.2a0.5O9VCn";
% Alpha_Mesh_Swc segmented mesh
swc_file = sprintf("swc_files/%s.swc",cellname); 
[results,femesh_cell,femesh_soma,femesh_neurites]= run_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
% Modified Ultraliser segmented mesh
[results_um,femesh_cell_um,femesh_soma_um,femesh_neurites_um]= run_simulations_neuron(mesh_um,setup_file,tetgen_options,swc_file);

signal_soma = real(results.mf_soma.signal/femesh_soma.total_volume);
signal_soma_um = real(results_um.mf_soma.signal/femesh_soma_um.total_volume);

bvals = results.setup.gradient.bvalues;
namplitude = length(bvals);

fig = figure;
subplot(2,1,1); hold on;
plot_neurites_soma(femesh_soma,femesh_neurites);hold on;
title('Alpha_Mesh_Swc','Interpreter','none');
subplot(2,1,2); hold on;
plot_neurites_soma(femesh_soma_um,femesh_neurites_um)
title('Modified Ultraliser');
sgtitle(sprintf("%s",cellname),'Interpreter','none');
 
r



