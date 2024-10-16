%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'));

setup_file='setup_pgse_microglia';
%% Get list of meshes
fid = fopen("cells_human_well_meshed.txt",'r');
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

for i =ncells
    mesh = meshes(i); type = types(i);
    [~,cellname,~] = fileparts(mesh); 
    mesh_um = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cellname);
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file =  sprintf("mesh_Files/Alzheimer_study/Soma/%s.ply",cellname);
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
                inds = {1,2,3,4,5,6,[7,8],9};
                inds_um = {1,[6],2,5,4,7,3,8};
                success = true;

            case "ctrl_010319_13_826-2_1"
                inds = {1,2,3,4,5,6,7,8};
                inds_um = {1,2,5,4,6,3,7,8};
                success = true;
                   
            otherwise
                try
                    p = align_dendrites(femesh_neurites,femesh_neurites_um);
                    inds = num2cell([1:length(femesh_neurites)]); inds_um = num2cell(p);
                    
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
                    plot_neurites_soma(femesh_soma_um,femesh_neurites_um)
                    title('Modified Ultraliser');
                    sgtitle(sprintf("%s",cellname),'Interpreter','none');
                end


        end
        if success
            plot_neurite_alignment(femesh_soma,femesh_neurites,inds,femesh_soma_um,femesh_neurites_um,inds_um,cellname);
            [process_signals,process_signals_um] = merge_process_signals(results,results_um,femesh_neurites,femesh_neurites_um,inds,inds_um);
            volumes_neurites = cellfun(@(i) sum(volumes_neurites(i)), inds);
            % volumes_neurites_um =  cellfun(@(i) sum(volumes_neurites_um(i)), inds_um);
            weighted_neurites = volumes_neurites'.*process_signals;
            weighted_neurites_um = volumes_neurites'.*process_signals;
            weighted_soma = real(results.mf_soma.signal);
            weighted_soma_um = signal_soma_um *femesh_soma.total_volume ;
            
            if length(inds) > 1
            cmpt_avg = (weighted_soma + sum(weighted_neurites,1))./femesh_cell.total_volume;
            cmpt_avg_um = (weighted_soma_um + sum(weighted_neurites_um,1))./femesh_cell.total_volume;
            else
                cmpt_avg = (weighted_soma' + weighted_neurites)./femesh_cell.total_volume;
                cmpt_avg_um = (weighted_soma_um' + weighted_neurites_um)./femesh_cell.total_volume;
            end
            
            save(sprintf('neuron_meshing_paper/microglia_output/%s_seg.mat',cellname),'signal_soma','signal_soma_um','process_signals','process_signals_um','bvals','cmpt_avg','cmpt_avg_um');     
        end
    end
end







%%

for i =11
    mesh = meshes(i); type = types(i);    
    [results_neuron,femesh,femesh_soma,femesh_dendrites]= load_mf_segmented_microglia(mesh,setup_file,tet,'1.0');
    [~,cellname,~] = fileparts(mesh); 
    mesh = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cellname);
    if isfile(mesh)
        icell = icell + 1;  meshes_tested(icell) = cellname;
        [results_neuron_ultraliser,femesh_ult,femesh_soma_ult,femesh_dendrites_ult] = load_mf_segmented_microglia(mesh,setup_file,tet,'1.0');
        bvals = results_neuron.setup.gradient.bvalues;
        % save(sprintf('neuron_meshing_paper/output_data/%s.mat',cell),'results_neuron','results_neuron_ultraliser','bvals');
        figdir="neuron_meshing_paper/aligning_microglia";
        fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
        fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
        savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cellname,results_neuron.setup.geometry.tetgen_options));
        savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cellname,results_neuron_ultraliser.setup.geometry.tetgen_options));
        signal = real(results_neuron.mf_soma.signal/results_neuron.soma_volume);
        signal_ultraliser = real(results_neuron_ultraliser.mf_soma.signal/results_neuron_ultraliser.soma_volume);
        n = length(bvals);
        error = abs((signal- signal_ultraliser));%./signal);
        soma_errors{icell} = error;
        % errors(icell) = error;
        disp(error);
        ndendrites = length(femesh_dendrites)
        signal_dend = zeros([ndendrites,size(signal)]);
        for iden = 1:ndendrites
            signal_dend(iden,:) = real(results_neuron.mf_dendrites{iden}.signal)./results_neuron.dendrite_volumes{iden}
        end
        ndendrites_ult = length(femesh_dendrites_ult)
        signal_dend_ultraliser = zeros([ndendrites,size(signal_ultraliser)]);
        for iden = 1:ndendrites_ult
            signal_dend_ultraliser(iden,:) = real(results_neuron_ultraliser.mf_dendrites{iden}.signal)./results_neuron_ultraliser.dendrite_volumes{iden};
        end
        if ndendrites == ndendrites_ult && cellname ~= "ctrl_150219_11_766-1_1" && cellname ~= "ctrl_08.03.19_5_714-5_1"
            disp("Can align");
            try
                p = align_dendrites(femesh_dendrites,femesh_dendrites_ult,true);
                if length(p) >1
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                end
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);
                
            catch 
                dendrite_errors{icell} = -1;
            end
        else
            dendrite_errors{icell} = -1;

            fprintf("Need to manually adjust cell %d\n",i);
            if cellname == "826_6_3"
                p = [5 3 1 7 2 4 6 8] %6 contains both 7 and 8 and 9, 10 are to be ignored.
                signal_dend(7,:) = (results_neuron.dendrite_volumes{7}*signal_dend(7,:) + ...
                    results_neuron.dendrite_volumes{8}*signal_dend(8,:))./(results_neuron.dendrite_volumes{8}+ results_neuron.dendrite_volumes{7});
                ind = [1 2 3 4 5 6 7 9];
                femesh_dendrites = femesh_dendrites(ind);
                signal_dend =  signal_dend(ind,:);results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(ind);
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);
            elseif cellname == "ctrl_150219_12_766-2_1"
                % Delete 12 14 from our mesh and also 10 11
                % Delete 9 11 12 from modified Ultraliser
                % 
                ind = [1 2 3 4 5 6 7 8 9 13];
                ind_ult = [1 2 3 4 5 6 7 8 10];
                femesh_dendrites = femesh_dendrites(ind); signal_dend = signal_dend(ind,:);
                results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(ind);
                femesh_dendrites_ult = femesh_dendrites_ult(ind_ult); 
                results_neuron_ultraliser.dendrite_volumes = results_neuron_ultraliser.dendrite_volumes(ind_ult);

                signal_dend_ultraliser = signal_dend_ultraliser(ind_ult,:);
                fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
                fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
                savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cellname,results_neuron.setup.geometry.tetgen_options));
                savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cellname,results_neuron_ultraliser.setup.geometry.tetgen_options));
                signal_dend(9,:)  = (femesh_dendrites{9}.total_volume*signal_dend(9,:) + ...
                    femesh_dendrites{10}.total_volume*signal_dend(10,:))./(femesh_dendrites{9}.total_volume + femesh_dendrites{10}.total_volume);
                signal_dend = signal_dend(1:9,:); femesh_dendrites = femesh_dendrites(1:9);                
                results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(1:9);

                p = [3 5 2 6 9 4 1 7 8];
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);

            elseif cellname ==  "ctrl_150219_14_766-4_1"
                % Remove dendrite 3.5 of ultraliser, 3 of our method
                ind = [1 2 4];
                ind_ult = [1 2 4];
                femesh_dendrites = femesh_dendrites(ind); signal_dend = signal_dend(ind,:);
                results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(ind);
                femesh_dendrites_ult = femesh_dendrites_ult(ind_ult); 
                results_neuron_ultraliser.dendrite_volumes = results_neuron_ultraliser.dendrite_volumes(ind_ult);

                signal_dend_ultraliser = signal_dend_ultraliser(ind_ult,:);
                fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
                fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
                savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cellname,results_neuron.setup.geometry.tetgen_options));
                savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cellname,results_neuron_ultraliser.setup.geometry.tetgen_options));
                p = [2 1 3];
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);
            elseif cellname ==    "ctrl_150219_15_766-5_1"
                % From my cell remove 1
                % From Ultraliser remove 6 , 9
                ind = 2:8;
                ind_ult = [1 2 3 4 5 7 8];
                femesh_dendrites = femesh_dendrites(ind); signal_dend = signal_dend(ind,:);
                results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(ind);
                femesh_dendrites_ult = femesh_dendrites_ult(ind_ult); 
                results_neuron_ultraliser.dendrite_volumes = results_neuron_ultraliser.dendrite_volumes(ind_ult);

                signal_dend_ultraliser = signal_dend_ultraliser(ind_ult,:);
                fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
                fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
                savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cellname,results_neuron.setup.geometry.tetgen_options));
                savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cellname,results_neuron_ultraliser.setup.geometry.tetgen_options));
                p = align_dendrites(femesh_dendrites,femesh_dendrites_ult,true)
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);

            elseif cellname == "ctrl_150219_11_766-1_1"
                ind = [1 2 3 4 5 6 8];
                ind_ult = [1 2 3 4 5 6 8];
                femesh_dendrites = femesh_dendrites(ind); signal_dend = signal_dend(ind,:);
                results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(ind);
                femesh_dendrites_ult = femesh_dendrites_ult(ind_ult); 
                results_neuron_ultraliser.dendrite_volumes = results_neuron_ultraliser.dendrite_volumes(ind_ult);

                signal_dend_ultraliser = signal_dend_ultraliser(ind_ult,:);
                fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
                fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
                savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cellname,results_neuron.setup.geometry.tetgen_options));
                savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cellname,results_neuron_ultraliser.setup.geometry.tetgen_options));
                p = align_dendrites(femesh_dendrites,femesh_dendrites_ult,true)
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);
            elseif cellname == "ctrl_08.03.19_5_714-5_1"
                ind = [1 2 4 5 6];
                ind_ult = [1 2 3 4 6];
                femesh_dendrites = femesh_dendrites(ind); signal_dend = signal_dend(ind,:);
                results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(ind);
                femesh_dendrites_ult = femesh_dendrites_ult(ind_ult); 
                results_neuron_ultraliser.dendrite_volumes = results_neuron_ultraliser.dendrite_volumes(ind_ult);

                signal_dend_ultraliser = signal_dend_ultraliser(ind_ult,:);
                fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
                fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
                savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cellname,results_neuron.setup.geometry.tetgen_options));
                savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cellname,results_neuron_ultraliser.setup.geometry.tetgen_options));
                p = align_dendrites(femesh_dendrites,femesh_dendrites_ult,true)
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);
            elseif cellname == "ctrl_010319_13_826-2_1"
                ind = 1:8;
                ind_ult = 1:8;
                femesh_dendrites = femesh_dendrites(ind); signal_dend = signal_dend(ind,:);
                results_neuron.dendrite_volumes = results_neuron.dendrite_volumes(ind);
                femesh_dendrites_ult = femesh_dendrites_ult(ind_ult); 
                results_neuron_ultraliser.dendrite_volumes = results_neuron_ultraliser.dendrite_volumes(ind_ult);
                signal_dend_ultraliser = signal_dend_ultraliser(ind_ult,:);
                fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
                fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
                savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cellname,results_neuron.setup.geometry.tetgen_options));
                savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cellname,results_neuron_ultraliser.setup.geometry.tetgen_options));
                p = align_dendrites(femesh_dendrites,femesh_dendrites_ult,true)
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser));%./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);

            end

            %
            % figdir="neuron_meshing_paper/aligning_microglia";
            % fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
            % fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
            % savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cell,results_neuron.setup.geometry.tetgen_options));
            % savefig(fig,sprintf('%s/%s_um_tet%s.fig',figdir,cell,results_neuron_ultraliser.setup.geometry.tetgen_options));
        end

    if cellname ~=   "ctrl_210219_10_781-5_1"
        results_neuron_ultraliser.dendrite_volumes = results_neuron_ultraliser.dendrite_volumes(p,:)
        total_volume = results_neuron.total_volume;
        soma_volume = results_neuron.soma_volume;
        dendrite_volumes = zeros(length(femesh_dendrites),1);
        total_volume_u = results_neuron_ultraliser.total_volume;
        soma_volume_u = results_neuron_ultraliser.soma_volume;
        dendrite_volumes_u = zeros(length(femesh_dendrites),1);
        
        for iden = 1:length(femesh_dendrites)
            dendrite_volumes(iden) =results_neuron.dendrite_volumes{iden};
            dendrite_volumes_u(iden) =results_neuron_ultraliser.dendrite_volumes{iden}

        end


        save(sprintf('neuron_meshing_paper/microglia_output/%s_seg.mat',cellname),'signal','signal_ultraliser','signal_dend','signal_dend_ultraliser','bvals','total_volume','total_volume_u','soma_volume','soma_volume_u','dendrite_volumes','dendrite_volumes_u');
    end
    end
end

%%
fid = fopen('output.txt','w');
for i = 1:length(meshes_tested)
    fprintf(fid,"%s\n",meshes_tested(i));
    fprintf(fid,"Max soma error = %f\n",max(soma_errors{i}));
    fprintf(fid,"Max dendrite error = %f\n",max(dendrite_errors{i}'));
    fprintf(fid,'\n');
end
