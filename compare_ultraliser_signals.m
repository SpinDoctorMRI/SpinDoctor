cd ..
addpath(genpath('src'))
addpath(genpath('setups'))
setup_file='setup_comparison';
tet="-pq1.2a0.5O9VCn";
tet_btpde="-pq1.2a0.1O9VCn";
% mesh='mesh_files/ultraliser/1-2-2-watertight.ply';
mesh='mesh_files/ultraliser_modified/1-2-2.CNG_um.ply';
results_neuron_ultraliser= load_mf_segmented(mesh,setup_file,tet,'3.0');
% btpde_ultraliser = load_btpde_cell(mesh,setup_file,tet_btpde);
% btpde_ultraliser

tet="-pq1.2a0.5O9VCn";
mesh='mesh_files/selected/1-2-2.CNG.ply';
results_neuron= load_mf_segmented(mesh,setup_file,tet,'3.0');
% btpde = load_btpde_cell(mesh,setup_file,tet_btpde);
% btpde

% relative_errors = abs(real(results_neuron.mf.signal/results_neuron.total_volume- results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume))./abs(real(results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume));
% relative_errors = abs(real(btpde.btpde.signal/btpde.total_volume- btpde_ultraliser.btpde.signal/btpde_ultraliser.total_volume))./abs(real(btpde_ultraliser.btpde.signal/btpde_ultraliser.total_volume));

% error = abs(real(btpde.btpde.signal- btpde_ultraliser.btpde.signal))./abs(real(btpde_ultraliser.btpde.signal));

bvals = results_neuron.setup.gradient.bvalues;

% relative_errors = squeeze(relative_errors);
set(groot,'defaultLineLineWidth',3.0)

% fig = figure(1);
% hold on;
% plot(bvals,relative_errors,'DisplayName','Errors BTPDE','Marker','o')
% grid on;xlabel('b-values');ylabel('Relative error');legend('Location','eastoutside','Interpreter','None');
% saveas(fig,'neuron_meshing_paper/ultraliser_modified_comparison.png')



% relative_errors = abs(real(results_neuron.mf.signal/results_neuron.total_volume- results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume))./abs(real(results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume));
bvals = results_neuron.setup.gradient.bvalues;

% relative_errors = squeeze(relative_errors);
set(groot,'defaultLineLineWidth',3.0)

% fig = figure(1);
% fig.Position(3:4) = [430, 266];
% fontsize(fig, 24, "points");set(gca,'FontName','cmr12');
% hold on;
% plot(bvals,relative_errors,'DisplayName','Errors MF','Marker','o')
% grid on;xlabel('b-values');ylabel('Relative error');legend('Location','eastoutside','Interpreter','None');
% saveas(fig,'neuron_meshing_paper/ultraliser_modified_comparison1','epsc')


% fig = figure(2);
% fig.Position(3:4) = [430, 266];
% fontsize(fig, 24, "points");set(gca,'FontName','cmr12');
% plot(bvals,results_neuron.mf.signal/results_neuron.total_volume,'DisplayName','Our method','Marker','o');
% hold on;
% plot(bvals,results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume,'DisplayName','Modified Ultraliser','Marker','x');
% grid on;xlabel('b-values');ylabel('Volume weighted signal');legend('Location','eastoutside','Interpreter','None');
% saveas(fig,'neuron_meshing_paper/ultraliser_modified_comparison2','epsc')

% fig = figure(3);
% fig.Position(3:4) = [430, 266];
% fontsize(fig, 24, "points");set(gca,'FontName','cmr12');
% plot(bvals,results_neuron.mf_soma.signal/results_neuron.soma_volume,'DisplayName','Our method','Marker','o');
% hold on;
% plot(bvals,results_neuron_ultraliser.mf_soma.signal/results_neuron_ultraliser.soma_volume,'DisplayName','Modified Ultraliser','Marker','x');
% grid on;xlabel('b-values');ylabel('Volume weighted signal');legend('Location','eastoutside','Interpreter','None');
% saveas(fig,'neuron_meshing_paper/ultraliser_modified_comparison_soma','png')

close all;

dendrite_mapping = {[3,5,7],[2,4],[6],[1],[9],[8]}
nseq = length(bvals);

fig = figure(7);
fig.Position(3:4) = [430, 266];
fontsize(fig, 24, "points");set(gca,'FontName','cmr12');

for i = 1:6
    fprintf('Dendrite %d\n',i);
    fig = figure(i);
    fig.Position(3:4) = [430, 266];
    fontsize(fig, 24, "points");set(gca,'FontName','cmr12');
    signals = real(results_neuron_ultraliser.mf_dendrites{i}.signal)/results_neuron_ultraliser.dendrite_volumes{i};
    ind = dendrite_mapping{i};
    relevent_results = [results_neuron.mf_dendrites{ind}];
    other_signals =  reshape(real([relevent_results.signal]),nseq,[]);
    other_signals = sum(other_signals,2);
    other_signals=other_signals./sum([results_neuron.dendrite_volumes{ind}]);
    plot(bvals,other_signals,'DisplayName','Our method','Marker','o');
    hold on;
    plot(bvals,signals,'DisplayName','Modified Ultraliser','Marker','x');
    grid on;xlabel('b-values');ylabel('Volume weighted signal');legend('Location','eastoutside','Interpreter','None');
    ylim([0,1]);
    saveas(fig,sprintf('neuron_meshing_paper/ultraliser_modified_comparison_dendrite_%d',i),'png')


    fig =figure(7);
    rel_errors = abs(signals-other_signals')./abs(signals)
    hold on;
    plot(bvals,rel_errors,'DisplayName',sprintf('Dendrite %d',i),'Marker','x');


end
grid on;xlabel('b-values');ylabel('Relative errors');
legend('Location','eastoutside','Interpreter','None');
saveas(fig,'neuron_meshing_paper/ultraliser_modified_comparison_dendrite_errors','png')

save('neuron_meshing_paper/output_data/comparison_data.mat','results_neuron','results_neuron_ultraliser')