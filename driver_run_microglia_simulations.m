%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))

setup_file='setup_pgse_microglia';
%% Get list of meshes
fid = fopen("cells_human.txt",'r');
meshes = textscan(fid,"%s");
fclose(fid);
meshes = meshes{1}; meshes = string(meshes);
file_parts = split(meshes,"/");
types = file_parts(:,3);

ncells = length(meshes);

meshes_um = replace(meshes,"Alzheimer_study","ultraliser_modified");
meshes_um = replace(meshes_um,".ply","_um.ply");
%% Run simulations for mesh_swc meshes
tetgen_options = "-pq1.2a0.5O9VCn";
for i = 1:ncells
    mesh = meshes(i);
    run_simulations_microglia(mesh,setup_file,tetgen_options);
end
%%
tetgen_options = "-pq1.2a0.05O9VCn";
for i = 1:ncells
    mesh = meshes(i);
    [~,cellname,~] = fileparts(mesh);
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    run_simulations_microglia(mesh,setup_file,tetgen_options,swc_file,soma_file);
end
%% Run simulations for modified ultraliser meshes
tetgen_options = "-pq1.2a0.5O9VCn";
for i = 1:ncells
    mesh = meshes_um(i); 
    run_simulations_microglia(mesh,setup_file,tetgen_options);
end
tetgen_options = "-pq1.2a0.05O9VCn";
for i = 1:ncells
    mesh = meshes_um(i);
    [~,cellname_um,~] = fileparts(mesh);
    cellname = replace(cellname_um,"_um","");
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    run_simulations_microglia(mesh,setup_file,tetgen_options,swc_file,soma_file);
end
%% Compare cell signals
tetgen_options = "-pq1.2a0.5O9VCn";

clear rel_errors abs_errors; icell = 0;
clear volumes; clear volumes_ult;
clear meshes_tested;
for i =1:ncells
    mesh = meshes(i); type = types(i);
    [~,cell,~] = fileparts(mesh); 
    mesh_um = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cell);
    if isfile(mesh_um)
        [results,femesh_cell,~,~]= load_simulations_microglia(mesh,setup_file,tetgen_options);
        [results_um,femesh_cell_um,~,~]= load_simulations_microglia(mesh_um,setup_file,tetgen_options);
        bvals = results.setup.gradient.bvalues;
        % Save direct comparison here for plots in
        % plotting_microglia.ipynb.
        signal = real(results.mf_cell.signal./femesh_cell.total_volume);
        signal_um = real(results_um.mf_cell.signal./femesh_cell_um.total_volume);

        signal = real(results.mf_cell.signal/femesh_cell.total_volume);
        signal_ultraliser = real(results_um.mf_cell.signal/femesh_cell_um.total_volume);
        namplitude = length(bvals);
        rel_error = abs((signal- signal_ultraliser)./signal);
        abs_error  = abs(signal- signal_ultraliser);
        icell = icell + 1; 
        rel_errors(icell,:) = rel_error; abs_errors(icell,:) = abs_error;
        meshes_tested(icell) = cell;
        volumes(icell) = femesh_cell.total_volume; volumes_ult(icell) = femesh_cell_um.total_volume;
        volume=femesh_cell.total_volume
        save(sprintf('neuron_meshing_paper/microglia_output/%s.mat',cell),'signal','signal_um','bvals','volume');

    end
end

fig = figure(1);
plot(bvals,rel_errors);
legend(meshes_tested,'Interpreter','none');
title("Relative errors");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

fig = figure(2);
plot(bvals,abs_errors);
legend(meshes_tested,'Interpreter','none');
title("Absolute errors");
grid on;
xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');


%% Compare segmented cell signals
disp("Running compare_microglia_segmented.m");
driver_compare_microglia_simulations_segmented;

