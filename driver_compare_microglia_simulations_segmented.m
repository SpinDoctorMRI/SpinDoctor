%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'));

setup_file='setup_pgse_microglia';
%% Get list of meshes
fid = fopen("cells_human.txt",'r');
meshes = textscan(fid,"%s");
fclose(fid);
meshes = meshes{1}; meshes = string(meshes);
file_parts = split(meshes,"/");
types = file_parts(:,3);

ncells = length(meshes);

%%
tet="-pq1.2a0.05O9VCn";


clear soma_errors; icell = 0;
clear dendrite_errors; clear meshes_tested;
% clear volumes; clear volumes_ult;
clear meshes_need_adjustment;

fid = fopen('neuron_meshing_paper/microglia_output/volumes.txt','w');
fprintf(fid,"\\hline \\\\ \\multirow{2}{*}{Cell} & \\multicolumn{2}{c|}{Alpha\\_ Mesh\\_ Swc} & \\multicolumn{2}{c|}{Modified Ultraliser} \\\\ \n")
fprintf(fid," & Soma fraction & Processes fraction & Soma fraction & Processes fraction \\\\ \n\\hline\n");

for i =1:ncells
    mesh = meshes(i); type = types(i);
    [~,cellname,~] = fileparts(mesh); 
    mesh_um = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cellname);
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file =  sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    if isfile(mesh_um)
        icell = icell + 1;
        meshes_tested(icell) = cellname;
        [results,femesh_cell,femesh_soma,femesh_neurites]= load_simulations_microglia(mesh,setup_file,tet,swc_file,soma_file);
        [results_um,femesh_cell_um,femesh_soma_um,femesh_neurites_um]= load_simulations_microglia(mesh_um,setup_file,tet,swc_file,soma_file);
        bvals = results.setup.gradient.bvalues;
        % Save direct comparison here for plots in
        % plotting_microglia.ipynb.
        signal_soma = real(results.mf_soma.signal/femesh_soma.total_volume);
        signal_soma_um = real(results_um.mf_soma.signal/femesh_soma_um.total_volume);
        namplitude = length(bvals);

        femesh_neurites_merged = [femesh_neurites{:}];
        femesh_neurites_merged_um = [femesh_neurites_um{:}];
        
        volumes_neurites = [femesh_neurites_merged.total_volume];
        volumes_neurites_um = [femesh_neurites_merged_um.total_volume];
        success = false;
        switch cellname

            case "826_6_3"
                % inds = {1,2,3,4,5,6,[7,8],9};
                % inds_um = {1,[6],2,5,4,7,3,8};
                inds = {1,[2,7],3,4,5,6,8,9};
                inds_um = {1,6,2,4,5,3,7,8};
                success = true;
            
            case "766_4_3"
                inds = {1,2,3,4,5};
                inds_um = {3,4,5,6,7};
                success = true;
            case "714_3_2"
                disp("Could remove first from both but not necessary");
                inds = num2cell(1:4);inds_um = num2cell(1:4); success = true;
            case "ctrl_010319_13_826-2_1"
                disp("Could remove first from both but not necessary");
                p = align_dendrites(femesh_neurites,femesh_neurites_um);
                inds = num2cell([1:length(femesh_neurites)]); inds_um = num2cell(p);
                
                success = true;
                   
            otherwise
                try
                    p = align_dendrites(femesh_neurites,femesh_neurites_um);
                    inds = num2cell([1:length(femesh_neurites)]); inds_um = num2cell(p);
                    
                    success = true;
              catch
                    meshes_need_adjustment(icell) = cellname;
                    fig = figure;
                    subplot(1,2,1); hold on;
                    plot_neurites_soma(femesh_soma,femesh_neurites);hold on;
                    title('Our method');
                    subplot(1,2,2); hold on;
                    plot_neurites_soma(femesh_soma_um,femesh_neurites_um)
                    title('Modified Ultraliser');
                    sgtitle(sprintf("%s",cellname),'Interpreter','none');
                    savefig(fig,sprintf('%s.fig',cellname));
                end


            end
        if success
            % plot_neurite_alignment(femesh_soma,femesh_neurites,inds,femesh_soma_um,femesh_neurites_um,inds_um,cellname);
            [process_signals,process_signals_um,volumes,volumes_um] = merge_process_signals(results,results_um,femesh_neurites,femesh_neurites_um,inds,inds_um);
            soma_volume = femesh_soma.total_volume;soma_volume_um = femesh_soma_um.total_volume;
            abs(process_signals-process_signals_um) 
            
            save(sprintf('neuron_meshing_paper/microglia_output/%s_seg.mat',cellname),'signal_soma','signal_soma_um','process_signals','process_signals_um','bvals','soma_volume','soma_volume_um','volumes','volumes_um');   
            fprintf(fid,"%s & %.2f & %.2f & %.2f & %.2f \\\\ \n",cellname,femesh_soma.total_volume/femesh_cell.total_volume,sum(volumes_neurites)/femesh_cell.total_volume,femesh_soma_um.total_volume/femesh_cell_um.total_volume,sum(volumes_neurites_um)/femesh_cell_um.total_volume)
            % fprintf(fid,"%s\n\n Our method\n Total volume %.2f\n Soma volume fraction %.2f\n Processes volume fraction %.2f\n",cellname,femesh_cell.total_volume,femesh_soma.total_volume/femesh_cell.total_volume,sum(volumes_neurites)/femesh_cell.total_volume);
            % fprintf(fid,"Modified Ultraliser\n Total volume %.2f\n Soma volume fraction %.2f\n Processes volume fraction %.2f\n",femesh_cell_um.total_volume,femesh_soma_um.total_volume/femesh_cell_um.total_volume,sum(volumes_neurites_um)/femesh_cell_um.total_volume);       
        else
            error()
        end
    end
end
fprintf(fid,'\\hline');
fclose(fid);
close all;