%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))

setup_file='setup_pgse_neuron_btpde';
setup_file_mf='setup_pgse_neuron';

%% Set mesh to examine

% mesh='mesh_files/selected/1-2-2.CNG.ply';
mesh='mesh_files/ultraliser_modified/1-2-2.CNG_um.ply';



%% Run simulations for mesh
tetgen_options_btpde = "-pq1.2a0.1O9VCn";
[~,cellname,~] = fileparts(mesh);
swc_file = sprintf("swc_files/%s.swc",cellname); 
run_simulations_neuron(mesh,setup_file,tetgen_options_btpde,swc_file);
run_simulations_neuron(mesh,setup_file,tetgen_options_btpde);



% %%
tetgen_options_mf = "-pq1.2a0.5O9VCn";
[~,cell,~] = fileparts(mesh); 
[results_btpde,femesh_cell_btpde,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options_btpde);
[results_mf,femesh_cell_mf,~,~]= load_simulations_neuron(mesh,setup_file_mf,tetgen_options_mf);
bvals = results_btpde.setup.gradient.bvalues;
% % Save direct comparison here for plots in
% % plotting_microglia.ipynb.
signal_btpde = real(results_btpde.btpde_cell.signal./femesh_cell_btpde.total_volume);
signal_mf = real(results_mf.mf_cell.signal./femesh_cell_mf.total_volume);

save(sprintf('neuron_meshing_paper/neuron_output/%s_convergence.mat',cell),'signal_btpde','signal_mf','bvals');
namplitude = length(bvals);
rel_error = abs((signal_btpde- signal_mf)./signal_btpde);
abs_error  = abs(signal_btpde- signal_mf);

fig = figure(1);
plot(bvals,rel_error);
title("Relative errors");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

fig = figure(2);
plot(bvals,abs_error);
title("Absolute errors");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

[results_btpde,femesh_cell_btpde,femesh_soma_btpde,femesh_neurites_btpde]= load_simulations_neuron(mesh,setup_file,tetgen_options_btpde,swc_file);
[results_mf,femesh_cell_mf,femesh_soma_mf,femesh_neurites_mf]= load_simulations_neuron(mesh,setup_file_mf,tetgen_options_mf,swc_file);

signal_soma_mf = real(results_mf.mf_soma.signal/femesh_soma_mf.total_volume);
signal_soma_btpde = real(results_btpde.btpde_soma.signal/femesh_soma_btpde.total_volume);
namplitude = length(bvals);
p = align_dendrites(femesh_neurites_mf,femesh_neurites_btpde);
inds = num2cell([1:length(femesh_neurites_mf)]); inds_btpde = num2cell(p);
fig = figure;
subplot(2,1,1); hold on;
plot_neurites_soma(femesh_soma_mf,femesh_neurites_mf);hold on;
title('h = 0.05');
subplot(2,1,2); hold on;
plot_neurites_soma(femesh_soma_btpde,femesh_neurites_btpde)
title('h = 0.01 ');
sgtitle(sprintf("%s",cellname),'Interpreter','none');
% savefig(fig,"neuron.fig");

plot_neurite_alignment(femesh_soma_mf,femesh_neurites_mf,inds,femesh_soma_btpde,femesh_neurites_btpde,inds_btpde,cellname);

femesh_neurites_btpde = femesh_neurites_btpde(p); btpde_neurites = results_btpde.btpde_neurites(p);
nneurites = length(femesh_neurites_btpde);
signal_neurites_btpde= zeros(nneurites,namplitude);
signal_neurites_mf= zeros(nneurites,namplitude);
for i = 1:nneurites
    signal_neurites_btpde(i,:) = real(btpde_neurites{i}.signal)/femesh_neurites_btpde{i}.total_volume;
    signal_neurites_mf(i,:) = real(results_mf.mf_neurites{i}.signal)/femesh_neurites_mf{i}.total_volume;
end
rel_error = abs(signal_soma_mf - signal_soma_btpde)./signal_soma_btpde;
abs_error = abs(signal_soma_mf - signal_soma_btpde);

neurites_rel_error=  abs((signal_neurites_btpde- signal_neurites_mf)./signal_neurites_btpde);
neurites_abs_error=  abs(signal_neurites_btpde- signal_neurites_mf);

% save(sprintf('neuron_meshing_paper/neuron_output/%s_convergence_seg.mat',cellname),'rel_error','abs_error','neurites_rel_error','neurites_abs_error','bvals');   
