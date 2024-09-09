function get_segmentation_info(mesh_path,tetgen_options)

    addpath(genpath('setups'));
    addpath(genpath('src'));

    % Launch setup here
    setup_comparison





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
    % exportgraphics(gcf,sprintf('neuron_meshing_paper/%s_tet%s.png',cellname,setup.geometry.tetgen_options),'BackgroundColor','none');
    savefig(fig,sprintf('neuron_meshing_paper/%s_tet%s.fig',cellname,setup.geometry.tetgen_options));
    total_volume = femesh.total_volume;
    soma_volume = femesh_soma.total_volume;
    dendrite_volume = 0;
    for iden = 1:length(femesh_dendrites)
        dendrite_volume = dendrite_volume + femesh_dendrites{iden}.total_volume;
    end
    fid = fopen('neuron_meshing_paper/output_data/segmentation_info.txt','a');
    fprintf(fid,'\\verb|%s| & %.1f & %.1f & %.1f \\\\ \n',cellname,total_volume,100*soma_volume/total_volume,100*dendrite_volume/total_volume);
    fclose(fid);
    quit