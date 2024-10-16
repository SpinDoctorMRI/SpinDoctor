function compare_ply(file1,file2) 
fprintf("Reading %s\n",file1);
[f,v,~] = read_ply(file1);
cell=triangulation(f,v);
fprintf("Reading %s\n",file2);
[f,v,~] = read_ply(file2);
cell2=triangulation(f,v);
[~,name1,~]=fileparts(file1);
[~,name2,~]=fileparts(file2);

fig=figure;
h=trisurf(cell,"FaceColor","b",'DisplayName',name1,"LineStyle","none");
set(h, "edgealpha", 0.5);
set(h, "facealpha", 0.5);
hold on;
h = trisurf(cell2,"FaceColor","r",'DisplayName',name2,"LineStyle","none");
set(h, "edgealpha", 0.5);
set(h, "facealpha", 0.5);
legend('Interpreter','none')
savefig(fig,sprintf('neuron_meshing_paper/output_data/%s_vs_%s.fig',name1,name2))