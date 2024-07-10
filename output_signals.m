function output_signals(meshname,setup_file,tetgen_options,ls,direc1,direc2)
% Launch setup here
run(sprintf("%s.m",setup_file));
fprintf("Running %s.m\n",setup_file)

setup.name = string(meshname);
setup.geometry.tetgen_options = tetgen_options;
setup.mf.length_scale = str2num(ls);
[setup, femesh, ~]  = prepare_simulation(setup);

setup.name=string(meshname);
[mesh_path,cellname,~] = fileparts(setup.name);

cell_root =split(cellname,{'_process','_soma'}); cell_root = cell_root{1};

save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);
save_path_cell = sprintf("%s/cell",save_path_root);
mf = load_mf(setup, save_path_cell, false);
if ~isdir(sprintf('%s/%s', direc1,cell_root))
    mkdir(sprintf('%s/%s', direc1,cell_root));
end
fid = fopen(sprintf('%s/%s/%s.signals',direc1,cell_root,cellname),'w');
fprintf(fid,'%f\n',real(mf.signal));
fclose(fid);

results = struct;
results.setup = setup;
results.mf = mf;
results.total_volume = femesh.total_volume;
if ~isdir(sprintf('%s', direc2))
    mkdir(sprintf('%s', direc2));
end
save(sprintf('%s/%s_cell.mat',direc2,cellname),'results');
disp('Completed')