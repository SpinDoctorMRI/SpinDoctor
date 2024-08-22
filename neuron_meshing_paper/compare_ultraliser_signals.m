cd ..
addpath(genpath('src'))
addpath(genpath('setups'))
setup_file='setup_bvalues_experiments'
tet="-pq1.1a1.0VCn";


mesh='mesh_files/ultraliser/1-2-1-watertight.ply'
results_neuron_ultraliser= load_mf_cell(mesh,setup_file,tet,'1.0')

% mesh='mesh_files/ultraliser/20-sn-1-watertight.ply'
% results_microglia_ultraliser= load_mf_cell(mesh,setup_file,tet,'1.0')

tet="-pq1.2a0.1VCn";
mesh='mesh_files/selected/1-2-1.CNG.ply'
results_neuron= load_mf_cell(mesh,setup_file,tet,'1.0')

% mesh='mesh_files/ultraliser/20-sn-1.CNG.ply'
% results_microglia= load_mf_cell(mesh,setup_file,tet,'1.0')

relative_errors = abs(real(results_neuron.mf.signal/results_neuron.total_volume- results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume))./abs(real(results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume));
bvals = results_neuron.setup.gradient.bvalues;

relative_errors = squeeze(relative_errors);
set(groot,'defaultLineLineWidth',3.0)

fig = figure(1)
plot(bvals,mean(relative_errors,2),'DisplayName','Direction-averaged error','Marker','x')
hold on;
plot(bvals,max(relative_errors'),'DisplayName','Direction-maximised error','Marker','o')
grid on;xlabel('b-values');ylabel('Relative error');legend('Location','eastoutside','Interpreter','None');
saveas(fig,'neuron_meshing_paper/ultraliser_comparison_1.png')

% save('neuron_meshing_paper/ultraliser_comparison_data.mat','relative_errors','bvals')
