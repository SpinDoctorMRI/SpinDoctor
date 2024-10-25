%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))

setup_file='setup_morez_microglia';
%% Get list of meshes
fid = fopen("cells_human.txt",'r');
meshes = textscan(fid,"%s");
fclose(fid);
meshes = meshes{1}; meshes = string(meshes);
file_parts = split(meshes,"/");
types = file_parts(:,3);

ncells = length(meshes);
%% Run simulations for mesh_swc meshes
% tetgen_options = "-pq1.2a0.5O9VCn";
% for i = 1:ncells
%     mesh = meshes(i);
%     run_simulations_microglia(mesh,setup_file,tetgen_options);
% end
% %%
% tetgen_options_seg = "-pq1.2a0.05O9VCn";
% for i = 1:ncells
%     mesh = meshes(i);
%     [~,cellname,~] = fileparts(mesh);
%     swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
%     run_simulations_microglia(mesh,setup_file,tetgen_options_seg,swc_file,soma_file);
% end
%% Compare cell signals
tetgen_options = "-pq1.2a0.5O9VCn";
tetgen_options_seg = "-pq1.2a0.05O9VCn";

for i =1%:ncells
    mesh = meshes(i); type = types(i);
    [~,cellname,~] = fileparts(mesh);
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    [results_cell,femesh_cell,~,~]= load_simulations_microglia(mesh,setup_file,tetgen_options);
    [results_seg,femesh_cell_seg,femesh_soma,femesh_processes] = load_simulations_microglia(mesh,setup_file,tetgen_options_seg,swc_file,soma_file);
    
    % Do any saving and processing here
end

