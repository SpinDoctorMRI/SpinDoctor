function results = driver_btpde(mesh_path,setup_file,tetgen_options,full_cell)
% driver_btpde  runs experiments.
% mesh_path : (str) path to mesh
% setup_file: (str) setup file to load
% tetgen_options : (str) tetgen_options for mesh. Defaults to '-pq1.2aVCn'
% ls: (float) length scale to determine eigenfunctions. Defaults to char_length_scale/5.
addpath(genpath('setups'));
addpath(genpath('src'));

% Launch setup here
run(sprintf("%s.m",setup_file));
fprintf("Running %s.m\n",setup_file)

% Flag f
full_cell = str2num(full_cell);



if nargin >= 3
    setup.geometry.tetgen_options=string(tetgen_options);
else
    fprintf('Applying default Tetgen parameters %s\n',setup.geometry.tetgen_options);
end

setup.name=string(mesh_path);
[mesh_path,cellname,~] = fileparts(setup.name);
[setup, femesh, ~, ~]  = prepare_simulation(setup);

save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);

if isfield(setup,'btpde')
    save_path_cell = sprintf("%s/cell",save_path_root);
    btpde_cell = solve_btpde(femesh, setup,save_path_cell,false);
    results.cell = btpde_cell;
    results.total_volume = femesh.total_volume;
    if full_cell ~= 1
        swc_file=sprintf("swc_files/%s.swc",cellname);
        tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
        [femesh_soma,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);
        % fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
        % saveas(fig,sprintf('figures/%s_tet%s.png',cellname,setup.geometry.tetgen_options));

        save_path_soma = sprintf("%s/soma",save_path_root);
        btpde_soma = solve_btpde(femesh_soma, setup,save_path_soma,false);
        ndendrites=length(femesh_dendrites);
        btpde_dendrites = cell(ndendrites,1);
        for i=1:ndendrites 
            save_path_dendrite = sprintf("%s/dendrite_%d",save_path_root,i);
            % Compute MF magnetization
            btpde_dendrites{i} = solve_btpde(femesh_dendrites{i}, setup,save_path_dendrite,false);
        end
        results.soma = btpde_soma;
        results.dendrites=btpde_dendrites;
        results.soma_volume = femesh_soma.total_volume;
        results.dendrite_volume = [femesh_dendrites{:}.total_volume];
    end


end
results.setup=setup;