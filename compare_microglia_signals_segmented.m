addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'))

setup_file='setup_pgse_ref_sol';

tet="-pq1.2a0.1O9VCn";
tet_btpde="-pq1.2a0.1O9VCn";

fid = fopen("cells_human.txt",'r');
tline = fgetl(fid);
i=1;
clear meshes;
while ischar(tline)
    meshes(i) = string(tline);
    tline = fgetl(fid);
    i = i + 1;
end
ncells = length(meshes)
clear types;

for i = 1:ncells
    file = meshes(i);
    sep_file = split(file,"/");
    types(i) = sep_file(3);
end
%%
clear soma_errors; icell = 0;
clear dendrite_errors; close all
clear meshes_tested;

% clear volumes; clear volumes_ult;
for i = 11
    mesh0 = meshes(i); type = types(i);    
    % [results_neuron,femesh,femesh_soma,femesh_dendrites]= load_mf_segmented_microglia(mesh,setup_file,tet,'1.0');

    [~,cell,~] = fileparts(mesh0); 
    cell
    mesh = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cell);
    if isfile(mesh)
        [results_neuron,femesh,femesh_soma,femesh_dendrites]= driver_mf_segmented_microglia(mesh0,setup_file,tet,'1.0');

        icell = icell + 1; meshes_tested(icell) = cell;
        % [results_neuron_ultraliser,femesh_ult,femesh_soma_ult,femesh_dendrites_ult] = load_mf_segmented_microglia(mesh,setup_file,tet,'1.0');
        [results_neuron_ultraliser,femesh_ult,femesh_soma_ult,femesh_dendrites_ult] = driver_mf_segmented_microglia(mesh,setup_file,tet,'1.0');

        bvals = results_neuron.setup.gradient.bvalues;
        % save(sprintf('neuron_meshing_paper/output_data/%s.mat',cell),'results_neuron','results_neuron_ultraliser','bvals');
        % figdir="neuron_meshing_paper/aligning_microglia";
        fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
        fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
        % savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cell,results_neuron.setup.geometry.tetgen_options));
        % savefig(fig_ult,sprintf('%s/%s_um_tet%s.fig',figdir,cell,results_neuron_ultraliser.setup.geometry.tetgen_options));
        signal = real(results_neuron.mf_soma.signal/results_neuron.soma_volume);
        signal_ultraliser = real(results_neuron_ultraliser.mf_soma.signal/results_neuron_ultraliser.soma_volume);
        n = length(bvals);
        error = abs((signal- signal_ultraliser)./signal);
        soma_errors{icell} = error;
        % errors(icell) = error;
        disp(error);
        ndendrites = length(femesh_dendrites)
        signal_dend = zeros([ndendrites,size(signal)]);
        for iden = 1:ndendrites
            signal_dend(iden,:) = real(results_neuron.mf_dendrites{iden}.signal)./results_neuron.dendrite_volumes(iden)
        end
        ndendrites_ult = length(femesh_dendrites_ult)
        signal_dend_ultraliser = zeros([ndendrites,size(signal_ultraliser)]);
        for iden = 1:ndendrites_ult
            signal_dend_ultraliser(iden,:) = real(results_neuron_ultraliser.mf_dendrites{iden}.signal)./results_neuron_ultraliser.dendrite_volumes(iden);
        end
        if i==2
            signal_dend_ultraliser = signal_dend_ultraliser(1:8,:);femesh_dendrites_ult = femesh_dendrites_ult(1:8);
            ndendrites_ult = 8;
        end
        if ndendrites == ndendrites_ult
            disp("Can align");
            try
                p = align_dendrites(femesh_dendrites(1:8),femesh_dendrites_ult)
                signal_dend_ultraliser = signal_dend_ultraliser(p,:);
                signal_dend = squeeze(signal_dend);signal_dend_ultraliser = squeeze(signal_dend_ultraliser);
                errors =  abs((signal_dend- signal_dend_ultraliser)./signal_dend);
                dendrite_errors{icell} = errors;
                disp(errors);
                
            end
        else
            fprintf("Need to manually adjust cell %d\n",i);
            % figdir="neuron_meshing_paper/aligning_microglia";
            % fig = plot_dendrites_soma(femesh_soma,femesh_dendrites);
            % fig_ult = plot_dendrites_soma(femesh_soma_ult,femesh_dendrites_ult);
            % savefig(fig,sprintf('%s/%s_tet%s.fig',figdir,cell,results_neuron.setup.geometry.tetgen_options));
            % savefig(fig,sprintf('%s/%s_um_tet%s.fig',figdir,cell,results_neuron_ultraliser.setup.geometry.tetgen_options));
        end
      plot_field_everywhere(femesh_soma,results_neuron.mf_soma.magnetization(10), 'My method')
    plot_field_everywhere(femesh_soma_ult,results_neuron_ultraliser.mf_soma.magnetization(10), 'Ultraliser')


    end
end
