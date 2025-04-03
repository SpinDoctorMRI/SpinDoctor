%% Add SpinDoctor to Path
restoredefaultpath;
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'))

setup_file='setup_pgse_neuron_final'; tetgen_options = "-pq1.2a0.5O9VCn";


% Default mesh level
mesh='mesh_files/selected/1-2-2.CNG.ply';
% Run full cell simulations
run_simulations_neuron(mesh,setup_file,tetgen_options);
Load and save full cell simulations
[results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
bvals = results.setup.gradient.bvalues;namplitude = length(bvals);
mf =results.mf_cell;
signals = real(mf.signal/femesh_cell.total_volume);
save("neuron_meshing_paper\neuron_output\1-2-2.CNG_signals.mat",'bvals','signals');

% Refinement level 1
mesh='mesh_files/selected/1-2-2.CNG_level1.ply';
% Run full cell simulations
run_simulations_neuron(mesh,setup_file,tetgen_options);
Load and save full cell simulations
[results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
bvals = results.setup.gradient.bvalues;namplitude = length(bvals);
mf =results.mf_cell;
signals = real(mf.signal/femesh_cell.total_volume);
save("neuron_meshing_paper\neuron_output\1-2-2.CNG_level1_signals.mat",'bvals','signals');

% Refinement level 2
mesh='mesh_files/selected/1-2-2.CNG_level2.ply';
% Run full cell simulations
run_simulations_neuron(mesh,setup_file,tetgen_options);
% Load and save full cell simulations
[results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
bvals = results.setup.gradient.bvalues;namplitude = length(bvals);
mf =results.mf_cell;
signals = real(mf.signal/femesh_cell.total_volume);
save("neuron_meshing_paper\neuron_output\1-2-2.CNG_level2_signals.mat",'bvals','signals');


%%Old section
% swc_file = "swc_files/1-2-2.CNG.swc";
% Run segmented cell simulations
% run_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);

%% Load and save segmented cell simulations
% [results_seg,femesh_cell_seg,femesh_soma,femesh_neurites]= load_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
% bvals = results_seg.setup.gradient.bvalues;namplitude = length(bvals);
% mf_soma =results_seg.mf_soma;
% signals_soma = real(mf_soma.signal/femesh_soma.total_volume);


% mf_neurites = [results_seg.mf_neurites{:}];
% nneurites = length(results_seg.mf_neurites);
% signals_neurites = real([mf_neurites.signal]);
% signals_neurites = reshape(signals_neurites,[namplitude,nneurites])';
% signal_neurite_compartment = sum(signals_neurites,1);
% volumes_neurites = [femesh_neurites.total_volume];
% volume_neurite_cmpt = sum(volumes_neurites);

% signals_neurites = signals_neurites/volume_neurite_cmpt;

% save("neuron_meshing_paper\neuron_output\1-2-2._seg.mat",'bvals','signals_soma','signal_neurites');

