%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))

setup_file='setup_pgse_neuron';
mesh='mesh_files/selected/1-2-2.CNG.ply';

%% Run simulations for mesh_swc mesh
% tetgen_options = "-pq1.2a0.5O9VCn";
tetgen_options = "-pq1.2a0.5O9VCn";
run_simulations_neuron(mesh,setup_file,tetgen_options);

[~,cellname,~] = fileparts(mesh);
swc_file = sprintf("swc_files/%s.swc",cellname); 
run_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
%% Run simulations for modified ultraliser meshes

run_simulations_neuron(mesh,setup_file,tetgen_options);

[~,cellname,~] = fileparts(mesh);
swc_file = sprintf("swc_files/%s.swc",cellname); 
run_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);

%% Compare cell signals
tetgen_options = "-pq1.2a0.5O9VCn";

[~,cell,~] = fileparts(mesh); 
mesh_um = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cell);
[results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
[results_um,femesh_cell_um,~,~]= load_simulations_neuron(mesh_um,setup_file,tetgen_options);
bvals = results.setup.gradient.bvalues;
% Save direct comparison here for plots in
% plotting_microglia.ipynb.
signal = real(results.mf_cell.signal./femesh_cell.total_volume);
signal_um = real(results_um.mf_cell.signal./femesh_cell_um.total_volume);

save(sprintf('neuron_meshing_paper/neuron_output/%s.mat',cell),'signal','signal_um','bvals');
signal = real(results.mf_cell.signal/femesh_cell.total_volume);
signal_ultraliser = real(results_um.mf_cell.signal/femesh_cell_um.total_volume);
namplitude = length(bvals);
rel_error = abs((signal- signal_ultraliser)./signal);
abs_error  = abs(signal- signal_ultraliser);

fig = figure(1);
plot(bvals,rel_error);
legend(meshes_tested,'Interpreter','none');
title("Relative errors");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

fig = figure(2);
plot(bvals,abs_error);
legend(meshes_tested,'Interpreter','none');
title("Absolute errors");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');


%% Compare segmented cell signals

[results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_neuron(mesh,setup_file,tetgen_options,swc_file);
[results_um,femesh_cell_um,femesh_soma_um,femesh_neurites_um]= load_simulations_neuron(mesh_um,setup_file,tetgen_options,swc_file);

% plotting_microglia.ipynb.
signal_soma = real(results.mf_soma.signal/femesh_soma.total_volume);
signal_soma_um = real(results_um.mf_soma.signal/femesh_soma_um.total_volume);
namplitude = length(bvals);

femesh_neurites_merged = [femesh_neurites{:}];
femesh_neurites_merged_um = [femesh_neurites_um{:}];

volumes_neurites = [femesh_neurites_merged.total_volume];
volumes_neurites_um = [femesh_neurites_merged_um.total_volume];
p = align_dendrites(femesh_neurites,femesh_neurites_um);
inds = num2cell([1:length(femesh_neurites)]); inds_um = num2cell(p);

            % meshes_need_adjustment(icell) = cellname;
            % figures =  findobj('type','figure');
            % nfigures = length(figures);
            % figure(1+nfigures);
            % subplot(1,2,1); hold on;
            % plot_neurites_soma(femesh_soma,femesh_neurites);hold on;
            % title('Our method');
            % subplot(1,2,2); hold on;
            % plot_neurites_soma(femesh_soma_um,femesh_neurites_um)
            % title('Modified Ultraliser');
            % sgtitle(sprintf("%s",cellname),'Interpreter','none');


plot_neurite_alignment(femesh_soma,femesh_neurites,inds,femesh_soma_um,femesh_neurites_um,inds_um,cellname);
[process_signals,process_signals_um] = merge_process_signals(results,results_um,femesh_neurites,femesh_neurites_um,inds,inds_um);
volumes_neurites = cellfun(@(i) sum(volumes_neurites(i)), inds);
% volumes_neurites_um =  cellfun(@(i) sum(volumes_neurites_um(i)), inds_um);
weighted_neurites = volumes_neurites'.*process_signals;
weighted_neurites_um = volumes_neurites'.*process_signals;
weighted_soma = real(results.mf_soma.signal);
weighted_soma_um = signal_soma_um *femesh_soma.total_volume ;

if length(inds) > 1
cmpt_avg = (weighted_soma + sum(weighted_neurites,1))./femesh_cell.total_volume;
cmpt_avg_um = (weighted_soma_um + sum(weighted_neurites_um,1))./femesh_cell.total_volume;
else
    cmpt_avg = (weighted_soma' + weighted_neurites)./femesh_cell.total_volume;
    cmpt_avg_um = (weighted_soma_um' + weighted_neurites_um)./femesh_cell.total_volume;
end

save(sprintf('neuron_meshing_paper/neuron_output/%s_seg.mat',cellname),'signal_soma','signal_soma_um','process_signals','process_signals_um','bvals','cmpt_avg','cmpt_avg_um');   

