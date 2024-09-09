function [results,femesh]= load_btpde_cell(mesh,setup_file,tetgen_options,reltol,abstol)
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
    if nargin >= 4 && ~isempty(str2num(reltol))
        setup.btpde.reltol = str2num(reltol);
        fprintf('Relative tolerance = %f\n',reltol);
    else
        fprintf('Applying default relative tolerance %f\n',setup.btpde.reltol);
    end

    if nargin == 5 && ~isempty(str2num(abstol))
        setup.btpde.abstol = str2num(abstol);
        fprintf('Absolute tolerance = %f\n',abstol);
    else
        fprintf('Applying default absolute tolerance %f\n',setup.btpde.abstol);
    end
    
    setup.name=string(mesh);
    [~,cellname,~] = fileparts(setup.name);
    [setup, femesh, ~, ~]  = prepare_simulation(setup);
    
    save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);
    
    if isfield(setup,'mf')
        save_path_cell = sprintf("%s/cell",save_path_root);
        %         % Compute MF magnetization
        btpde = load_btpde(setup,save_path_cell,false);
    end
    results.setup = setup; results.btpde = btpde;
    results.total_volume = femesh.total_volume;
    results.mesh_size=length(femesh.points{1});
    toc