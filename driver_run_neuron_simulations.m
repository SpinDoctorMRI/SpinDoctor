%% Add SpinDoctor to Path
restoredefaultpath;
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'))

setup_file='setup_pgse_neuron'; tetgen_options = "-pq1.2a0.5O9VCn";
mesh='mesh_files/selected/1-2-2.CNG.ply';
[~,cellname,~] = fileparts(mesh);

mesh_um = sprintf("mesh_files/ultraliser_modified/%s_um.ply",cellname);

%% Run segmented cell simulations
% Alpha_Mesh_Swc segmented mesh

swc_file = sprintf("swc_files/%s.swc",cellname); 
[results,femesh_cell,femesh_soma,femesh_neurites]= run_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
% Modified Ultraliser segmented mesh
[results_um,femesh_cell_um,femesh_soma_um,femesh_neurites_um]= run_simulations_neuron(mesh_um,setup_file,tetgen_options,swc_file);

%% Load segmented cell simulations
swc_file = sprintf("swc_files/%s.swc",cellname); 
[results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
[results_um,femesh_cell_um,femesh_soma_um,femesh_neurites_um]= load_simulations_neuron(mesh_um,setup_file,tetgen_options,swc_file);


%% Plot Meshes
fig = figure;
subplot(2,1,1); hold on;
plot_neurites_soma(femesh_soma,femesh_neurites);hold on;
title('Alpha_Mesh_Swc','Interpreter','none');
subplot(2,1,2); hold on;
plot_neurites_soma(femesh_soma_um,femesh_neurites_um)
title('Modified Ultraliser');
sgtitle(sprintf("%s",cellname),'Interpreter','none');

%%
% Soma signal and volumes
signal_soma = real(results.mf_soma.signal);
signal_soma_um = real(results_um.mf_soma.signal);
volume_soma = femesh_soma.total_volume;
volume_soma_um = femesh_soma_um.total_volume;

% Plot soma signals
figure;plot(bvals,signal_soma/volume_soma); hold on;
plot(bvals,signal_soma_um/volume_soma_um);
title("Comparing volume weighted soma signals");
ylabel("Volume weighted signals");
legend(["Alpha_Mesh_Swc","Modified Ultraliser"],'Interpreter','none');
grid on;xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');



% %% Access neurite signals and volumes

bvals = results.setup.gradient.bvalues;namplitude = length(bvals);


% Alpha_Mesh_Swc
mf_neurites = [results.mf_neurites{:}];
nneurites = length(results.mf_neurites);
signals_neurites = real([mf_neurites.signal]);
signals_neurites = reshape(signals_neurites,[namplitude,nneurites])';
signal_neurite_compartment = sum(signals_neurites,1);

% Modified Ultraliser
mf_neurites_um = [results_um.mf_neurites{:}];
nneurites_um = length(results_um.mf_neurites);
signals_neurites_um = real([mf_neurites_um.signal]);
signals_neurites_um = reshape(signals_neurites_um,[namplitude,nneurites_um])';
signal_neurite_compartment_um = sum(signals_neurites_um,1);

% Volumes neurites
volume = femesh_cell.total_volume;
volume_um = femesh_cell_um.total_volume;
femesh_n = [femesh_neurites{:}];
volumes_neurites = [femesh_n.total_volume];
volume_neurite_cmpt = sum(volumes_neurites);
femesh_n = [femesh_neurites_um{:}];
volumes_neurites_um = [femesh_n.total_volume];
volume_neurite_cmpt_um = sum(volumes_neurites_um);

clear femesh_n;

% Plot neurite compartments
figure;
plot(bvals,signal_neurite_compartment/volume); hold on;
plot(bvals,signal_neurite_compartment_um/volume_um);
title("Comparing volume weighted neurite compartment signals");
ylabel("Volume weighted signals");
legend(["Alpha_Mesh_Swc","Modified Ultraliser"],'Interpreter','none');
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

% Plot cell from weighted compartments

signal_cell = (signal_neurite_compartment + signal_soma)/volume;
signal_cell_um = (signal_neurite_compartment_um + signal_soma_um)/volume_um;

figure;
plot(bvals,signal_cell); hold on;
plot(bvals,signal_cell_um);
title("Comparing volume weighted signals");
ylabel("Volume weighted signals");
legend(["Alpha_Mesh_Swc","Modified Ultraliser"],'Interpreter','none');
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

figure;
plot(bvals,signal_cell,'LineStyle','-','Color','b','LineWidth',2); hold on;
plot(bvals,signal_cell_um,'LineStyle','--','Color','b','LineWidth',2);
plot(bvals,signal_soma/volume_soma,'LineStyle','-','Color','g','LineWidth',2);
plot(bvals,signal_soma_um/volume_soma_um,'LineStyle','--','Color','g','LineWidth',2);
plot(bvals,signal_neurite_compartment/volume_neurite_cmpt,'LineStyle','-','Color','magenta','LineWidth',2);
plot(bvals,signal_neurite_compartment_um/volume_neurite_cmpt_um,'LineStyle','--','Color','magenta','LineWidth',2);
title("Comparing volume weighted signals");
ylabel("Volume weighted signals");
legend(["Alpha_Mesh_Swc: cell","Modified Ultraliser: cell", ...
    "Alpha_Mesh_Swc: soma","Modified Ultraliser: soma", ...
    "Alpha_Mesh_Swc: neurites","Modified Ultraliser: neurites"], ...
    'Interpreter','none');
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

figure;
plot(bvals,abs(signal_cell-signal_cell_um),'LineStyle','-','Color','b','LineWidth',2); hold on;
plot(bvals,abs(signal_neurite_compartment/volume_neurite_cmpt-signal_neurite_compartment_um/volume_neurite_cmpt_um), ...
    'LineStyle','-','Color','magenta','LineWidth',2);
plot(bvals,abs(signal_soma/volume_soma-signal_soma_um/volume_soma_um),'LineStyle','-','Color','g','LineWidth',2);
title("Absolute difference of volume weighted signals");
ylabel("Volume weighted signals");
legend(["Cell","Soma","Neurites"], ...
    'Interpreter','none');
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');




%% Save signals and volumes
% fprintf("Alpha_Mesh_Swc \n Modified Ultraliser volume = %f.\n",volume,volume_um);
fprintf("Alpha_Mesh_Swc:\nCell volume = %f, Soma volume = %f, Neurites volume = %f.\n",volume,volume_soma,volume_neurite_cmpt);
fprintf("Modified Ultraliser:\nCell volume = %f, Soma volume = %f, Neurites volume = %f.\n",volume_um,volume_soma_um,volume_neurite_cmpt_um);


save(sprintf('neuron_meshing_paper/neuron_output/%s_seg.mat',cellname), ...
    'signal_soma','signal_soma_um','signals_neurites','signals_neurites_um', ...
    'bvals','volume_soma','volume_soma_um','volumes_neurites', ...
    'volumes_neurites_um',"signal_cell","signal_cell_um");   



%% Run simulations for full cell

% % Alpha_Mesh_Swc mesh
% [results,femesh_cell,~,~]=run_simulations_neuron(mesh,setup_file,tetgen_options);
% 
% % Modified Ultraliser mesh
% [results_um,femesh_cell_um,~,~]= run_simulations_neuron(mesh_um,setup_file,tetgen_options);
% 
% %% Load simulations
% 
% % Alpha_Mesh_Swc mesh
% [results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
% 
% % Modified Ultraliser mesh
% [results_um,femesh_cell_um,~,~]= load_simulations_neuron(mesh_um,setup_file,tetgen_options);
% 
% %% Visualise and compare signals
% 
% % Extract bvalues and signals
% bvals = results.setup.gradient.bvalues;
% signal = real(results.mf_cell.signal./femesh_cell.total_volume);
% signal_um = real(results_um.mf_cell.signal./femesh_cell_um.total_volume);
% 
% 
% volume=femesh_cell.total_volume;
% volume_um = femesh_cell_um.total_volume;
% 
% filename=sprintf('neuron_meshing_paper/neuron_output/%s.mat',cellname);
% fprintf('Saving signals and volume to %s\n',filename);
% save(filename,'signal','signal_um','bvals','volume','volume_um');
% namplitude = length(bvals);
% 
% rel_diff = abs((signal- signal_um)./signal);
% abs_diff= abs(signal- signal_um);
% 
% figure;
% 
% figure;plot(bvals,signal);hold on;plot(bvals,signal_um);
% title("Comparing volume weighted signals");
% ylabel("Volume weighted signals");
% legend(["Alpha_Mesh_Swc","Modified Ultraliser"],'Interpreter','none');
% grid on;
% xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');
% 
% 
% figure;plot(bvals,rel_diff);title("Absolute relative difference");grid on;
% xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');
% 
% figure;plot(bvals,abs_diff);title("Absolute difference");grid on;
% xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');
% 
% fprintf("Alpha_Mesh_Swc volume = %f,\n Modified Ultraliser volume = %f.\n",volume,volume_um);
% fprintf("Here we only plot and consider the volume-weighted signals:\n")
% fprintf("Maximum relative signal difference = %f,\n Maximum absolute signal difference =%f.\n",max(rel_diff),max(abs_diff));
% 


