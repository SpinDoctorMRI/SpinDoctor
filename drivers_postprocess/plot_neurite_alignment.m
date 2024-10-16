function plot_neurite_alignment(femesh_soma,femesh_neurites,inds,femesh_soma_um,femesh_neurites_um,inds_um,cellname)
ngroups = length(inds_um);
figures =  findobj('type','figure');
nfigures = length(figures);
for igroup = 1:ngroups
    fig = figure('WindowState', 'maximized');
    subplot(1,2,1); hold on;
    plot_neurites_soma(femesh_soma,femesh_neurites,inds{igroup});hold on;
    title('Our method');
    subplot(1,2,2); hold on;
    plot_neurites_soma(femesh_soma_um,femesh_neurites_um,inds_um{igroup})
    title('Modified Ultraliser');
    sgtitle(sprintf("%s, process group %d",cellname,igroup),'Interpreter','none');

end