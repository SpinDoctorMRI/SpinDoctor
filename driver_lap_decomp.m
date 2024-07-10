function results = driver_lap_decomp(mesh_path,exp_str,tetgen_options,ls)
% driver  runs experiments.
% mesh_path : (str) path to mesh
% exp_str: (str) descriptor for experiments. Copy of results will be saved to signals/cell_name_tetgen_options_exp_str
% tetgen_options : (str) tetgen_options for mesh. Defaults to '-pq1.2aVCn'
% ls: (float) length scale to determine eigenfunctions. Defaults to char_length_scale/5.
if nargin == 1
    exp_str=''
end

addpath(genpath('setups'));
addpath(genpath('src'));

%Change setup here
setup_morez;




if nargin >= 3
    setup.geometry.tetgen_options=string(tetgen_options);
else
    fprintf('Applying default Tetgen parameters %s\n',setup.geometry.tetgen_options);
end
if nargin == 4
    setup.mf.length_scale = str2num(ls);
    fprintf('Eigenvalue length scale = %f',setup.mf.length_scale);
else
    fprintf('Applying default Eigenvalue length scale %s\n',setup.geometry.tetgen_options);
end

setup.name=string(mesh_path);
[mesh_path,cellname,~] = fileparts(setup.name);

disp('Segmenting cell geometry')
swc_file=sprintf("swc_files/%s.swc",cellname);
[setup, femesh, ~, ~]  = prepare_simulation(setup);
tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
[femesh_soma,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);
% fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
% saveas(fig,sprintf('figures/%s_tet%s.png',cellname,setup.geometry.tetgen_options));

save_path_root=sprintf('saved_simul/%s_tet%s',cellname,setup.geometry.tetgen_options);


disp("Begin laplacian decomposition of components");

if isfield(setup,'mf')
    disp("Decomposing cell")
    save_path_cell = sprintf("%s/cell",save_path_root);
    try
    compute_laplace_eig(femesh, setup.pde, setup.mf,save_path_cell);  
    catch
        fid = fopen('problem_lap_eigs.txt','a');
        fprintf(fid,'%s ls = %f\n',save_path_cell,setup.mf.length_scale);
        fclose(fid)
    end
    disp("Decomposing soma")       
    save_path_soma = sprintf("%s/soma",save_path_root);
    try
    compute_laplace_eig(femesh_soma, setup.pde, setup.mf,save_path_soma);         
    catch
        fid = fopen('problem_lap_eigs.txt','a');
        fprintf(fid,'%s ls = %f\n',save_path_soma,setup.mf.length_scale);
        fclose(fid)
    end
    ndendrites=length(femesh_dendrites);
    for i=1:ndendrites 
        fprintf('Decomposing dendrite %d/%d\n',i,ndendrites);
        save_path_dendrite = sprintf("%s/dendrite_%d",save_path_root,i);
        try
        compute_laplace_eig(femesh_dendrites{i}, setup.pde, setup.mf,save_path_dendrite);      
        catch
            fid = fopen('problem_lap_eigs.txt','a');
            fprintf(fid,'%s ls = %f\n',save_path_dendrite,setup.mf.length_scale);
            fclose(fid);
        end   
    end
end
disp('Complete');