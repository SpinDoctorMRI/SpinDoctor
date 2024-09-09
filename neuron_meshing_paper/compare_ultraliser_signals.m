cd ..
addpath(genpath('src'))
addpath(genpath('setups'))
setup_file='setup_comparison';
tet="-pq1.2a0.5O9VCn";
set(groot,'defaultLineLineWidth',3.0)

%%
mesh='mesh_files/ultraliser_modified/1-2-2.CNG_um.ply';
[mf_ultraliser,femesh_u]= load_mf_cell(mesh,setup_file,tet,'3.0');


%%
mesh='mesh_files/selected/1-2-2.CNG.ply';
[mf_neuron,femesh]= load_mf_cell(mesh,setup_file,tet,'3.0');
% [btpde_neuron,femesh_b]= load_btpde_cell(mesh,setup_file,tet_btpde);
%%
relative_errors = abs(real(mf_neuron.mf.signal/mf_neuron.total_volume- mf_ultraliser.mf.signal/mf_ultraliser.total_volume))./abs(real(mf_ultraliser.mf.signal/mf_ultraliser.total_volume));

bvals = mf_neuron.setup.gradient.bvalues;
relative_errors = squeeze(relative_errors);
%%

fig = figure(1);
plot(bvals,relative_errors,'DisplayName','Matrix formalism relative signal differences','Marker','x')
grid on;xlabel('b-values');ylabel('Relative error');legend('Location','eastoutside','Interpreter','None');


fig = figure(2);
plot(bvals,mf_neuron.mf.signal/mf_neuron.total_volume,'DisplayName','Our method','Marker','x')
hold on;
plot(bvals,mf_ultraliser.mf.signal/mf_ultraliser.total_volume,'DisplayName','Modified Ultraliser','Marker','x')

grid on;xlabel('b-values');ylabel('Volume weighted signal');legend('Location','eastoutside','Interpreter','None');



%%
set(0,'defaulttextinterpreter','latex')

ifield = 5;
plot_field_everywhere(femesh, mf_neuron.mf.magnetization(ifield), sprintf('Our method MF b=%f',bvals(ifield)))
set(gca,'fontsize',24)

plot_field_everywhere(femesh_u, mf_ultraliser.mf.magnetization(ifield), sprintf('Modified Ultraliser MF b=%f',bvals(ifield)))
set(gca,'fontsize',24)













