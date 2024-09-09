function [results,femesh]= load_mf_segmented(mesh,setup_file,tetgen_options,ls)
    % driver_cell  runs experiments for only the full cell.
    % mesh_path : (str) path to mesh
    % exp_str: (str) descriptor for experiments. Copy of results will be saved to signals/cell_name_tetgen_options_exp_str
    % tetgen_options : (str) tetgen_options for mesh. Defaults to '-pq1.2aVCn'
    % ls: (float) length scale to determine eigenfunctions. Defaults to char_length_scale/5.
    tic
    addpath(genpath('setups'));
    addpath(genpath('src'));
    
    % Launch setup here
    run(sprintf("%s.m",setup_file));
    fprintf("Running %s.m\n",setup_file)
    
    
    if nargin >=3
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
    
    setup.name=string(mesh);
    [mesh_path,cellname,~] = fileparts(setup.name);
    swc_file=sprintf("swc_files/%s.swc",cellname);

    [setup, femesh, ~, ~]  = prepare_simulation(setup);
    tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
    [femesh_soma,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);
    save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);
    fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
    savefig(fig,sprintf('neuron_meshing_paper/figures/%s_tet%s.fig',cellname,setup.geometry.tetgen_options));

    if isfield(setup,'mf')
        save_path_cell = sprintf("%s/cell",save_path_root);
        %         % Compute MF magnetization
        mf = load_mf(setup,save_path_cell,false);
        save_path_soma = sprintf("%s/soma",save_path_root);
        mf_soma = load_mf(setup,save_path_soma,false);
        ndendrites=length(femesh_dendrites);
        mf_dendrites = cell(ndendrites,1);
        for i=1:ndendrites 
            save_path_dendrite = sprintf("%s/dendrite_%d",save_path_root,i);
            mf_dendrites{i} = load_mf(setup,save_path_dendrite,false);
        end
    end
    results.setup = setup; results.mf = mf;
    results.mf_soma = mf_soma;
    results.mf_dendrites = mf_dendrites;
    results.total_volume = femesh.total_volume;
    results.soma_volume = femesh_soma.total_volume;
    results.dendrite_volumes = cell(ndendrites,1);
    for i = 1:ndendrites
        results.dendrite_volumes{i} =femesh_dendrites{i}.total_volume;
    end
    results.mesh_size=length(femesh.points{1});
    toc