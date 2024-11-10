%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'))

setup_file='setup_pgse_neuron';
mesh='mesh_files/selected/1-2-2.CNG.ply';

%% Run simulations for mesh_swc mesh
tetgen_options = "-pq1.2a0.5O9VCn";
% run_simulations_neuron(mesh,setup_file,tetgen_options);
% [~,cellname,~] = fileparts(mesh);
% swc_file = sprintf("swc_files/%s.swc",cellname); 
% run_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
% %% Run simulations for modified ultraliser meshes
tetgen_options = "-pq1.2a0.5O9VCn";
[~,cellname,~] = fileparts(mesh);

mesh_um = sprintf("mesh_files/ultraliser_modified/%s_um.ply",cellname);
% run_simulations_neuron(mesh_um,setup_file,tetgen_options);

swc_file = sprintf("swc_files/%s.swc",cellname); 
% run_simulations_neuron(mesh_um,setup_file,tetgen_options,swc_file);

%% Compare cell signals

% [~,cell,~] = fileparts(mesh); 
% mesh_um = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cell);
[results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
[results_um,femesh_cell_um,~,~]= load_simulations_neuron(mesh_um,setup_file,tetgen_options);
bvals = results.setup.gradient.bvalues;
% % Save direct comparison here for plots in
% % plotting_microglia.ipynb.
signal = real(results.mf_cell.signal./femesh_cell.total_volume);
signal_um = real(results_um.mf_cell.signal./femesh_cell_um.total_volume);

volume=femesh_cell.total_volume;
save(sprintf('neuron_meshing_paper/neuron_output/%s.mat',cellname),'signal','signal_um','bvals','volume');namplitude = length(bvals);
rel_diff = abs((signal- signal_um)./signal);
abs_diff= abs(signal- signal_um);

fig = figure(1);
plot(bvals,rel_diff);
% legend(meshes_tested,'Interpreter','none');
title("Relative absolute difference");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

fig = figure(2);
plot(bvals,abs_diff);
% legend(meshes_tested,'Interpreter','none');
title("Absolute difference");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');


%% Compare segmented cell signals

[results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
[results_um,femesh_cell_um,femesh_soma_um,femesh_neurites_um]= load_simulations_neuron(mesh_um,setup_file,tetgen_options,swc_file);

% % plotting_microglia.ipynb.
signal_soma = real(results.mf_soma.signal/femesh_soma.total_volume);
signal_soma_um = real(results_um.mf_soma.signal/femesh_soma_um.total_volume);
namplitude = length(bvals);

% femesh_neurites_merged = [femesh_neurites{:}];
% femesh_neurites_merged_um = [femesh_neurites_um{:}];

% volumes_neurites = [femesh_neurites_merged.total_volume];`
% volumes_neurites_um = [femesh_neurites_merged_um.total_volume];
inds = {1,[2,4],[3,5,7],9,8,6};
inds_um = {4,2,1,5,6,3};
% p = align_dendrites(femesh_neurites,femesh_neurites_um);
% inds = num2cell([1:length(femesh_neurites)]); inds_um = num2cell(p);

fig = figure;
subplot(2,1,1); hold on;
plot_neurites_soma(femesh_soma,femesh_neurites);hold on;
title('Our method');
subplot(2,1,2); hold on;
plot_neurites_soma(femesh_soma_um,femesh_neurites_um)
title('Modified Ultraliser');
sgtitle(sprintf("%s",cellname),'Interpreter','none');
savefig(fig,"neuron.fig");


% plot_neurite_alignment(femesh_soma,femesh_neurites,inds,femesh_soma_um,femesh_neurites_um,inds_um,cellname);
[process_signals,process_signals_um,volumes,volumes_um] = merge_process_signals(results,results_um,femesh_neurites,femesh_neurites_um,inds,inds_um);
soma_volume = femesh_soma.total_volume;soma_volume_um = femesh_soma_um.total_volume;
signal_differnce = abs(process_signals-process_signals_um) ;
save(sprintf('neuron_meshing_paper/neuron_output/%s_seg.mat',cellname),'signal_soma','signal_soma_um','process_signals','process_signals_um','bvals','soma_volume','soma_volume_um','volumes','volumes_um');   


% save(sprintf('neuron_meshing_paper/neuron_output/%s_seg.mat',cellname),'signal_soma','signal_soma_um','process_signals','process_signals_um','bvals');   

