function results= load_mf_cell(mesh,setup_file,tetgen_options,ls)
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
    [~,cellname,~] = fileparts(setup.name);
    [setup, femesh, ~, ~]  = prepare_simulation(setup);
    
    save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);
    
    if isfield(setup,'mf')
        save_path_cell = sprintf("%s/cell",save_path_root);
        %         % Compute MF magnetization
        mf = load_mf(setup,save_path_cell,false);
    end
    results.setup = setup; results.mf = mf;
    results.total_volume = femesh.total_volume;
    results.mesh_size=length(femesh.points{1});
    toc