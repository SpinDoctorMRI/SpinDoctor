function results = driver_segmented(mesh_path,setup_file,tetgen_options,ls)
% driver_segmented  runs experiments for all components of cell.
% mesh_path : (str) path to mesh
% exp_str: (str) descriptor for experiments. Copy of results will be saved to signals/cell_name_tetgen_options_exp_str
% tetgen_options : (str) tetgen_options for mesh. Defaults to '-pq1.2aVCn'
% ls: (float) length scale to determine eigenfunctions. Defaults set in setup.

addpath(genpath('setups'));
addpath(genpath('src'));

% Launch setup here
run(sprintf("%s.m",setup_file));
fprintf("Running %s.m\n",setup_file)





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


% Create finite element meshes
setup.name=string(mesh_path);
[mesh_path,cellname,~] = fileparts(setup.name);
swc_file=sprintf("swc_files/%s.swc",cellname);
[setup, femesh, ~, ~]  = prepare_simulation(setup);
tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
[femesh_soma,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);
fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
saveas(fig,sprintf('figures/%s_tet%s.png',cellname,setup.geometry.tetgen_options));

save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);

if isfield(setup,'mf')
    % Simulations on full cell
    save_path_cell = sprintf("%s/cell",save_path_root);
    lap_eig = compute_laplace_eig(femesh, setup.pde, setup.mf,save_path_cell);         
    % Compute MF magnetization
    mf = solve_mf(femesh, setup, lap_eig,save_path_cell,false);
    
    % Simulations on soma
    save_path_soma = sprintf("%s/soma",save_path_root);
    lap_eig = compute_laplace_eig(femesh_soma, setup.pde, setup.mf,save_path_soma);         
    % Compute MF magnetization
    mf_soma = solve_mf(femesh_soma, setup, lap_eig,save_path_soma,false);

    %Simulations on dendrites
    ndendrites=length(femesh_dendrites);
    mf_dendrites = cell(ndendrites,1);
    for i=1:ndendrites 
        save_path_dendrite = sprintf("%s/dendrite_%d",save_path_root,i);
        lap_eig = compute_laplace_eig(femesh_dendrites{i}, setup.pde, setup.mf,save_path_dendrite);         
        % Compute MF magnetization
        mf_dendrites{i} = solve_mf(femesh_dendrites{i}, setup, lap_eig,save_path_dendrite,false);
    end
    % Store results
    results.cell = mf;
    results.soma = mf_soma;
    results.dendrites=mf_dendrites;
    results.total_volume = femesh.total_volume;
    results.soma_volume = femesh_soma.total_volume;
    results.dendrite_volume = [femesh_dendrites{:}.total_volume];
end
results.setup=setup;