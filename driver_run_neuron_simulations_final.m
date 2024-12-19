%% Add SpinDoctor to Path
restoredefaultpath;
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'))

setup_file='setup_pgse_neuron_final'; tetgen_options = "-pq1.2a0.5O9VCn";
mesh='mesh_files/selected/1-2-2.CNG.ply';
[~,cellname,~] = fileparts(mesh);

%% Run segmented cell simulations
% Alpha_Mesh_Swc segmented mesh

swc_file = sprintf("swc_files/%s.swc",cellname); 
[results,femesh_cell,~,~]= run_simulations_neuron(mesh,setup_file,tetgen_options);

%% Load segmented cell simulations
swc_file = sprintf("swc_files/%s.swc",cellname); 
[results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
bvals = results.setup.gradient.bvalues;namplitude = length(bvals);
%%
mf =results.mf_cell;

%% Plot Meshes
fig = figure;
plot_field_everywhere(femesh_cell,mf.magnetization,"");
%%
signals = real(mf.signal/femesh_cell.total_volume);
save("neuron_meshing_paper\neuron_output\1-2-2.CNG_signals.mat",'bvals','signals');