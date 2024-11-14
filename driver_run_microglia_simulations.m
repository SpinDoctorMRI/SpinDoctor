%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))
setup_file='setup_pgse_microglia';

fid = fopen("cells_human.txt",'r');
meshes = textscan(fid,"%s");
fclose(fid);
meshes = meshes{1}; meshes = string(meshes);
file_parts = split(meshes,"/");
types = file_parts(:,3);

ncells = length(meshes);

meshes_um = replace(meshes,"Alzheimer_study","ultraliser_modified");
meshes_um = replace(meshes_um,".ply","_um.ply");
%% Run simulations for Alpha_Mesh_Swc meshes
tetgen_options = "-pq1.2a0.05O9VCn";
for i = 1:ncells
    mesh = meshes(i);
    [~,cellname,~] = fileparts(mesh);
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    run_simulations_microglia(mesh,setup_file,tetgen_options,swc_file,soma_file);
end
%% Run simulations for Modified Ultraliser meshes
tetgen_options = "-pq1.2a0.05O9VCn";
for i = 1:ncells
    mesh = meshes_um(i);
    [~,cellname_um,~] = fileparts(mesh);
    cellname = replace(cellname_um,"_um","");
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    run_simulations_microglia(mesh,setup_file,tetgen_options,swc_file,soma_file);
end


%% Compare segmented cell signals

tet="-pq1.2a0.05O9VCn";

% Control y axis of plots
ymax = 0.08;




fid = fopen('neuron_meshing_paper/microglia_output/volumes.txt','w');
% fprintf(fid,"\\hline \\\\ \\multirow{2}{*}{Cell} & \\multicolumn{3}{c|}{{\\it Alpha\\_Mesh\\_Swc}}" + ...
    % " & \\multicolumn{3}{c}{{\\it Modified Ultraliser}} \\\\ \n");
fprintf(fid,"\\hline \n \\multirow{2}{*}{\\textbf{Cell}} & \\multirow{2}{*}{\\textbf{Method}} &" + ...
    "\\multicolumn{2}{c|}{\\textbf{Volume}} &\\multirowcell{2}{\\textbf{Soma volume}\\\\ \\textbf{fraction}} \\\\ \n " + ...
    " & & \\textbf{Soma} &\\multicolumn{1}{c|}{\\textbf{Processes}} & \\\\ \\hline \n");
    % fprintf(fid," & Soma & Processes & Soma fraction & Soma & Processes & Soma fraction \\\\ \n\\hline\n");

for i =1:ncells
    mesh = meshes(i); type = types(i);
    [~,cellname,~] = fileparts(mesh); 
    mesh_um = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cellname);
    swc_file = sprintf("swc_files/%s.swc",cellname); 
    soma_file =  sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    [results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_microglia(mesh,setup_file,tet,swc_file,soma_file);
    [results_um,femesh_cell_um,femesh_soma_um,femesh_neurites_um]= load_simulations_microglia(mesh_um,setup_file,tet,swc_file,soma_file);
    bvals = results.setup.gradient.bvalues; namplitude = length(bvals);
    
    % Soma signal and volumes
    signal_soma = real(results.mf_soma.signal);
    signal_soma_um = real(results_um.mf_soma.signal);
    volume_soma = femesh_soma.total_volume;
    volume_soma_um = femesh_soma_um.total_volume;

    % Alpha_Mesh_Swc Neurites
    mf_neurites = [results.mf_neurites{:}];
    nneurites = length(results.mf_neurites);
    signals_neurites = real([mf_neurites.signal]);
    signals_neurites = reshape(signals_neurites,[namplitude,nneurites])';
    signal_neurite_compartment = sum(signals_neurites,1);
    
    % Modified Ultraliser Neurites
    mf_neurites_um = [results_um.mf_neurites{:}];
    nneurites_um = length(results_um.mf_neurites);
    signals_neurites_um = real([mf_neurites_um.signal]);
    signals_neurites_um = reshape(signals_neurites_um,[namplitude,nneurites_um])';
    signal_neurite_compartment_um = sum(signals_neurites_um,1);
    
    % Volumes neurites
    volume = femesh_cell.total_volume;
    volume_um = femesh_cell_um.total_volume;
    femesh_n = [femesh_neurites{:}];
    volumes_neurites = [femesh_n.total_volume];
    volume_neurite_cmpt = sum(volumes_neurites);
    femesh_n = [femesh_neurites_um{:}];
    volumes_neurites_um = [femesh_n.total_volume];
    volume_neurite_cmpt_um = sum(volumes_neurites_um);
    clear femesh_n;

    % Plot cell from weighted compartments
    
    signal_cell = (signal_neurite_compartment + signal_soma)/volume;
    signal_cell_um = (signal_neurite_compartment_um + signal_soma_um)/volume_um;

    figure;
    plot(bvals,abs(signal_cell-signal_cell_um),'LineStyle','-','Color','b','LineWidth',2); hold on;
    plot(bvals,abs(signal_neurite_compartment/volume_neurite_cmpt-signal_neurite_compartment_um/volume_neurite_cmpt_um), ...
        'LineStyle','-','Color','magenta','LineWidth',2);
    plot(bvals,abs(signal_soma/volume_soma-signal_soma_um/volume_soma_um),'LineStyle','-','Color','g','LineWidth',2);
    title(sprintf("Absolute difference of volume weighted signals: %s",cellname),'Interpreter','none');
    ylabel("Volume weighted signals");
    legend(["Cell","Soma","Neurites"], ...
        'Interpreter','none');
    grid on;
    xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');
    ylim([0,ymax]);

    save(sprintf('neuron_meshing_paper/microglia_output/%s_seg.mat',cellname),'signal_cell','signal_cell_um','signal_soma','signal_soma_um', ...
        'signal_neurite_compartment','signal_neurite_compartment_um','bvals','volume_soma','volume_soma_um', ...
        'volume_neurite_cmpt','volume_neurite_cmpt_um','volume','volume_um');   
    % fprintf(fid,"%s & %.2g & %.2g & %.2g & %.2g & %.2g & %.2g \\\\ \n",cellname, ...
    %     volume_soma,volume_neurite_cmpt, volume_soma/volume ,...
    %     volume_soma_um,volume_neurite_cmpt_um, volume_soma_um/volume_um);
    fprintf(fid,"\\multirow{2}{*}{\\verb|%s|} & {\\it AMS} & %.2g & %.2g & %.2f \\\\ \n" + ...
        "& {\\it MU} & %.2g & %.2g & %.2f \\\\ \n \\hline",cellname, ...
        volume_soma,volume_neurite_cmpt, volume_soma/volume ,...
        volume_soma_um,volume_neurite_cmpt_um, volume_soma_um/volume_um);
end
fclose(fid);



%% Run simulations for full cell
% tetgen_options = "-pq1.2a0.5O9VCn";
% 
% % Alpha_Mesh_Swc
% for i = 1:ncells
%     mesh = meshes(i);
%     run_simulations_microglia(mesh,setup_file,tetgen_options);
% end
% 
% % Modified Ultraliser
% for i = 1:ncells
%     mesh = meshes_um(i); 
%     run_simulations_microglia(mesh,setup_file,tetgen_options);
% end
%% Compare signals for full cell
% tetgen_options = "-pq1.2a0.5O9VCn";
% 
% clear rel_errors abs_errors; icell = 0;
% clear volumes; clear volumes_ult;
% clear meshes_tested;
% for i =1:ncells
%     mesh = meshes(i); type = types(i);
%     [~,cell,~] = fileparts(mesh); 
%     mesh_um = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cell);
%     if isfile(mesh_um)
%         [results,femesh_cell,~,~]= load_simulations_microglia(mesh,setup_file,tetgen_options);
%         [results_um,femesh_cell_um,~,~]= load_simulations_microglia(mesh_um,setup_file,tetgen_options);
%         bvals = results.setup.gradient.bvalues;
%         % Save direct comparison here for plots in
%         % plotting_microglia.ipynb.
%         signal = real(results.mf_cell.signal./femesh_cell.total_volume);
%         signal_um = real(results_um.mf_cell.signal./femesh_cell_um.total_volume);
% 
%         signal = real(results.mf_cell.signal/femesh_cell.total_volume);
%         signal_ultraliser = real(results_um.mf_cell.signal/femesh_cell_um.total_volume);
%         namplitude = length(bvals);
%         rel_error = abs((signal- signal_ultraliser)./signal);
%         abs_error  = abs(signal- signal_ultraliser);
%         icell = icell + 1; 
%         rel_errors(icell,:) = rel_error; abs_errors(icell,:) = abs_error;
%         meshes_tested(icell) = cell;
%         volumes(icell) = femesh_cell.total_volume; volumes_ult(icell) = femesh_cell_um.total_volume;
%         volume=femesh_cell.total_volume;
%         save(sprintf('neuron_meshing_paper/microglia_output/%s.mat',cell),'signal','signal_um','bvals','volume');
% 
%     end
% end
% 
% figure;plot(bvals,rel_errors);
% legend(meshes_tested,'Interpreter','none');
% title("Relative absolute difference of volume weighted signals");
% grid on;
% xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');
% 
% figure;plot(bvals,abs_errors);
% legend(meshes_tested,'Interpreter','none');
% title("Absolute difference of volume weighted signals");
% grid on;
% xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

