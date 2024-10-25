%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'));

setup_file_btpde='setup_pgse_microglia_btpde';
setup_file_mf='setup_pgse_microglia';
%% Get list of meshes
fid = fopen("cells_human.txt",'r');
meshes = textscan(fid,"%s");
fclose(fid);
meshes = meshes{1}; meshes = string(meshes);
file_parts = split(meshes,"/");
types = file_parts(:,3);

ncells = length(meshes);

%%
tetgen_options_btpde="-pq1.1a0.01O9VCn";
tetgen_options_mf="-pq1.2a0.05O9VCn";


clear soma_errors; icell = 0;
clear dendrite_errors; clear meshes_tested;
% clear volumes; clear volumes_ult;
clear meshes_need_adjustment;


for i =1%:ncells
    mesh = meshes(i); type = types(i);
    [~,cellname,~] = fileparts(mesh); 
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file =  sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    icell = icell + 1;
    meshes_tested(icell) = cellname;
    % [results_mf,femesh_cell_mf,femesh_soma_mf,femesh_neurites_mf]= load_simulations_microglia(mesh,setup_file_mf,tetgen_options_mf,swc_file,soma_file);
    % [results_btpde,femesh_cell_btpde,femesh_soma_btpde,femesh_neurites_btpde]= load_simulations_microglia(mesh,setup_file_btpde,tetgen_options_btpde,swc_file,soma_file);
    bvals = results_mf.setup.gradient.bvalues;
    % Save direct comparison here for plots in
    % plotting_microglia.ipynb.
    signal_soma_mf = real(results_mf.mf_soma.signal/femesh_soma_mf.total_volume);
    signal_soma_btpde = real(results_btpde.btpde_soma.signal/femesh_soma_btpde.total_volume);
    namplitude = length(bvals);

    femesh_neurites_merged_mf = [femesh_neurites_mf{:}];
    femesh_neurites_merged_btpde = [femesh_neurites_btpde{:}];
    
    volumes_neurites_mf = [femesh_neurites_merged_mf.total_volume];
    volumes_neurites_btpde = [femesh_neurites_merged_btpde.total_volume];
    success = false;
    switch cellname

        % case "826_6_3"
        %     inds = {1,2,3,4,5,6,[7,8],9};
        %     inds_btpde = {1,[6],2,5,4,7,3,8};
        %     success = true;

        % case "ctrl_010319_13_826-2_1"
        %     inds = {1,2,3,4,5,6,7,8};
        %     inds_btpde = {1,2,5,4,6,3,7,8};
        %     success = true;
                
        otherwise
            try
                p = align_dendrites(femesh_neurites_mf,femesh_neurites_btpde);
                inds_mf = num2cell([1:length(femesh_neurites)]); inds_btpde = num2cell(p);
                
                success = true;
            catch
                meshes_need_adjustment(icell) = cellname;
                figures =  findobj('type','figure');
                nfigures = length(figures);
                figure(1+nfigures);
                subplot(1,2,1); hold on;
                plot_neurites_soma(femesh_soma,femesh_neurites);hold on;
                title('Our method');
                subplot(1,2,2); hold on;
                plot_neurites_soma(femesh_soma_btpde,femesh_neurites_btpde)
                title('Modified Ultraliser');
                sgtitle(sprintf("%s",cellname),'Interpreter','none');
            end


    end
    if success
        % plot_neurite_alignment(femesh_soma,femesh_neurites,inds,femesh_soma_btpde,femesh_neurites_btpde,inds_btpde,cellname);
        % [process_signals,process_signals_btpde] = merge_process_signals(results_mf,results_btpde,femesh_neurites_mf,femesh_neurites_btpde,inds_mf,inds_btpde);
        
        % save(sprintf('neuron_meshing_paper/microglia_output/%s_seg.mat',cellname),'signal_soma','signal_soma_btpde','process_signals','process_signals_btpde','bvals','cmpt_avg','cmpt_avg_btpde');       
    end
end
% fclose(fid);
close all;