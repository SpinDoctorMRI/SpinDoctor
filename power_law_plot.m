function power_law_plot(meshname,direc,tetgen_options,ls) 
cd ..
addpath(genpath("src"))
addpath(genpath("setups"))
set(0, 'DefaultLineLineWidth', 2);

setup_bvalues_experiments;

if nargin >= 3
    setup.geometry.tetgen_options=string(tetgen_options);
else
    fprintf('Applying default Tetgen parameters %s\n',setup.geometry.tetgen_options);
end
if nargin == 4 && ~isempty(str2num(ls))
    setup.mf.length_scale = str2num(ls);
    fprintf('Eigenvalue length scale = %f\n',setup.mf.length_scale);
else
    fprintf('Applying default Eigenvalue length scale %f\n',setup.mf.length_scale);
end

setup.name=string(meshname);
[mesh_path,cellname,~] = fileparts(setup.name);
swc_file=sprintf("swc_files/%s.swc",cellname);
[setup, femesh, ~, ~]  = prepare_simulation(setup);
tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
[femesh_soma,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);

save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);
save_path_cell = sprintf("%s/cell",save_path_root);
mf_cell = load_mf(setup,save_path_cell,false);
save_path_soma = sprintf("%s/cell",save_path_root);
mf_soma = load_mf(setup,save_path_soma,false);
ndendrites = length(femesh_dendrites);
mf_dends = cell(ndendrites,1);
for i = 1:ndendrites
    save_path_dendrite = sprintf("%s/dendrite_%d",save_path_root,i);
    mf_dends{i}=  load_mf(setup,save_path_dendrite,false);
end

fig = figure(1);
x = 1./sqrt(setup.gradient.values(1:end-2));
direction_avged_cell = real(mean(mf_cell.signal/femesh.total_volume,4));
direction_avged_cell = squeeze(direction_avged_cell);
scatter(x,direction_avged_cell(1:end-2),60,'o','DisplayName','Cell signal','LineWidth',2.5);
hold on
a = polyfit(x,direction_avged_cell(1:end-2),1);
plot(x,a(1)*x + a(2),'--','DisplayName',sprintf('S = %f b^{-1/2} + %f',a(1),a(2)));
direction_avged_cell = real(mean(mf_soma.signal/femesh_soma.total_volume,4));
direction_avged_cell = squeeze(direction_avged_cell);
scatter(x,direction_avged_cell(1:end-2),60,'o','DisplayName','Soma signal','LineWidth',2.5);
hold on
a = polyfit(x,direction_avged_cell(1:end-2),1);
plot(x,a(1)*x + a(2),'--','DisplayName',sprintf('S = %f b^{-1/2} + %f',a(1),a(2)));




ndendrites = length(mf_dends);
cmap = colormap(jet);
stepsize = floor(size(cmap,1)/ndendrites);
colors = cmap((1:ndendrites)*stepsize,:);

for iden = 1:ndendrites
    x = 1./sqrt(setup.gradient.values(1:end-2));
    direction_avged_cell = real(mean(mf_dends{iden}.signal/femesh_dendrites{iden}.total_volume,4));
    direction_avged_cell = squeeze(direction_avged_cell);
    scatter(x,direction_avged_cell(1:end-2),60,colors(iden,:),'x','DisplayName',sprintf('dendrite %d',iden),'LineWidth',2.5);
    hold on
    a = polyfit(x,direction_avged_cell(1:end-2),1);
    plot(x,a(1)*x + a(2),'Color',colors(iden,:),'LineStyle','--','DisplayName',sprintf('S = %f b^{-1/2} + %f',a(1),a(2)));
end
xlabel('b^{-1/2} [mm \cdot s^{-1/2}]')
ylabel('direction-averaged signal')
% title(sprintf('%s',cellname),'Interpreter','none')
legend;
legend('Location','northeastoutside')
% ylim([0,1])
a = gca;
a.FontSize=12; 
a.LineWidth = 1.5
grid on;
saveas(fig,sprintf('%s/%s_power_law.png',direc,cellname));
fprintf('Saved to %s/%s_power_law.png\n',direc,cellname);
quit