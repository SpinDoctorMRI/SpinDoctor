%% Add SpinDoctor to Path
addpath(genpath('src'))
addpath(genpath('setups'))
addpath(genpath('drivers_postprocess'))

setup_file_btpde='setup_pgse_microglia_btpde';
setup_file_mf='setup_pgse_microglia';

% Get list of meshes
fid = fopen("cells_human_well_meshed.txt",'r');
meshes = textscan(fid,"%s");
fclose(fid);
meshes = meshes{1}; meshes = string(meshes);
file_parts = split(meshes,"/");
types = file_parts(:,3);

ncells = length(meshes);

meshes_um = replace(meshes,"Alzheimer_study","ultraliser_modified");
meshes_um = replace(meshes_um,".ply","_um.ply");
%% Run simulations for mesh_swc meshes
tetgen_options="-pq1.1a0.01O9VCn";
for i = 1%:ncells
    mesh = meshes(i);
    % run_simulations_microglia(mesh,setup_file_btpde,tetgen_options);
    [~,cellname,~] = fileparts(mesh);
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    % run_simulations_microglia(mesh,setup_file_btpde,tetgen_options,swc_file,soma_file);
end
%% Run simulations for modified ultraliser meshes
tetgen_options = "-pq1.10.05O9VCn";
% tetgen_options = "-pq1.2a0.05O9VCn";
for i = 1:ncells
    mesh = meshes_um(i); 
    % run_simulations_microglia(mesh,setup_file_btpde,tetgen_options);
    [~,cellname_um,~] = fileparts(mesh);
    cellname = replace(cellname_um,"_um","");
    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    % run_simulations_microglia(mesh,setup_file_btpde,tetgen_options,swc_file,soma_file);
end
%% Compare cell signals
tetgen_options_mf = "-pq1.2a0.5O9VCn";
tetgen_options_btpde = "-pq1.1a0.01O9VCn";
clear rel_errors abs_errors; icell = 0;
clear volumes; clear volumes_ult;
clear meshes_tested;
clear cellnames ;
% for i =1:ncells
    % mesh = meshes(i); type = types(i);
    % [~,cellname,~] = fileparts(mesh); 
    % mesh = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cellname);
    % [results_btpde,femesh_cell_btpde,~,~]= load_simulations_microglia(mesh,setup_file_btpde,tetgen_options_btpde);
    % [results_mf,femesh_cell_mf,~,~]= load_simulations_microglia(mesh,setup_file_mf,tetgen_options_mf);
    % bvals = results_mf.setup.gradient.bvalues;
    % % Save direct comparison here for plots in
    % % plotting_microglia.ipynb.
    % signal_btpde = real(results_btpde.btpde_cell.signal./femesh_cell_btpde.total_volume);
    % signal_mf = real(results_mf.mf_cell.signal./femesh_cell_mf.total_volume);
    % % signal_btpde = real(results_btpde.btpde_cell.signal);
    % % signal_mf = real(results_mf.mf_cell.signal); 
    % namplitude = length(bvals);
    % rel_error = abs((signal_btpde- signal_mf)./signal_btpde);
    % abs_error  = abs(signal_btpde- signal_mf);
    % save(sprintf('neuron_meshing_paper/microglia_output/%s_convergence.mat',cellname),'rel_error','abs_error','bvals');

    % icell = icell + 1; 
    % rel_errors(icell,:) = rel_error; abs_errors(icell,:) = abs_error;
    % cellnames(i) = cellname;
% end

% fig = figure(1);
% plot(bvals,rel_errors);
% legend(cellnames,'Interpreter','none');
% title("Relative errors");
% grid on;
% xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');

% fig = figure(2);
% plot(bvals,abs_errors);
% legend(cellnames,'Interpreter','none');
% title("Absolute errors");
% grid on;
% xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');


%% Compare segmented cell signals
close all;
tetgen_options_mf = "-pq1.2a0.05O9VCn";
tetgen_options_btpde = "-pq1.1a0.01O9VCn";
clear rel_errors abs_errors; icell = 0;
clear volumes; clear volumes_ult;
clear meshes_tested;
clear cellnames ;

clear meshes_need_adjustment; imistake=0;
for i =1%:ncells
    mesh = meshes(i); type = types(i);
    [~,cellname,~] = fileparts(mesh); 
    % mesh = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cellname);

    swc_file = sprintf("swc_files/%s.swc",cellname); soma_file = sprintf("mesh_files/Alzheimer_study/Soma/%s.ply",cellname);
    [results_btpde,femesh_cell_btpde,femesh_soma_btpde,femesh_neurites_btpde]= load_simulations_microglia(mesh,setup_file_btpde,tetgen_options_btpde,swc_file,soma_file);
    [results_mf,femesh_cell_mf,femesh_soma_mf,femesh_neurites_mf]= load_simulations_microglia(mesh,setup_file_mf,tetgen_options_mf,swc_file,soma_file);
    bvals = results_mf.setup.gradient.bvalues;
    soma_btpde = real(results_btpde.btpde_soma.signal)./femesh_soma_btpde.total_volume;
    soma_mf = real(results_mf.mf_soma.signal)./femesh_soma_mf.total_volume; 
    success = false;
    [~,cellname,~] = fileparts(mesh); 

    switch cellname
        case "826_6_3"
            inds_btpde = {1,2,3,4,5,6,7,8};
            inds_mf = {1,3,4,5,6,[7,2],8,9};
            success =true;
            btpde_neurites_data =[results_btpde.btpde_neurites{:}];
            btpde_neurites = zeros(nneurites,namplitude);
            mf_neurites = zeros(nneurites,namplitude);
            if nneurites > 1
                for i = 1:nneurites
                    btpde_neurites(i,:) = btpde_neurites_data(i).signal/femesh_neurites_btpde{i}.total_volume;
                    
                end
            else
                btpde_neurites = btpde_neurites_data.signal/femesh_neurites_btpde{1}.total_volume;          
            end

            mf_neurites_data = [results_mf.mf_neurites{:}];
            for ii = 1:nneurites
                S = zeros(1,namplitude);
                V = 0
                for jj = 1: length(inds_mf{ii})
                    j = inds_mf{ii}(jj);
                    V = V + femesh_neurites_mf{j}.total_volume;
                    S = S + mf_neurites_data(j).signal;
                end
                mf_neurites(ii,:) = real(S/V);
            end

        otherwise
            try
                p = align_dendrites(femesh_neurites_btpde,femesh_neurites_mf);
                inds_btpde = num2cell([1:length(femesh_neurites_btpde)]); inds_mf = num2cell(p);
                success= true;
                btpde_neurites_data =[results_btpde.btpde_neurites{:}];
                mf_neurites_data = [results_mf.mf_neurites{p}];
                btpde_neurites = zeros(nneurites,namplitude);
                mf_neurites = zeros(nneurites,namplitude);
                
                femesh_neurites_mf = femesh_neurites_mf(p);
                if nneurites > 1
                    for i = 1:nneurites
                        mf_neurites(i,:) = mf_neurites_data(i).signal/femesh_neurites_mf{i}.total_volume;
                        btpde_neurites(i,:) = btpde_neurites_data(i).signal/femesh_neurites_btpde{i}.total_volume;
                        
                    end
                else
                    mf_neurites = mf_neurites_data.signal/femesh_neurites_mf{1}.total_volume;
                    btpde_neurites = btpde_neurites_data.signal/femesh_neurites_btpde{1}.total_volume;          
                end
            catch
                imistake = imistake + 1;
                meshes_need_adjustment(imistake) = cellname;
                fig =figure;
                subplot(1,2,1); hold on;
                plot_neurites_soma(femesh_soma_btpde,femesh_neurites_btpde);hold on;
                title('BTPDE');
                subplot(1,2,2); hold on;
                plot_neurites_soma(femesh_soma_mf,femesh_neurites_mf)
                title('MF');
                sgtitle(sprintf("%s",cellname),'Interpreter','none');
                savefig(fig,sprintf('%s.fig',cellname));
            end
        end
    if success
        % plot_neurite_alignment(femesh_soma_btpde,femesh_neurites_btpde,inds_btpde,femesh_soma_mf,femesh_neurites_mf,inds_mf,cellname);
        namplitude=length(bvals); nneurites = length(femesh_neurites_btpde);
       
        
        % btpde_neurites =reshape(real([btpde_neurites.signal]),[nneurites,namplitude]);
        % mf_neurites =reshape(real([mf_neurites.signal]),[nneurites,namplitude]);
        % btpde_neurites = zeros(nneurites,namplitude);
        % mf_neurites = zeros(nneurites,namplitude);
        
        % femesh_neurites_mf = femesh_neurites_mf(p);
        % if nneurites > 1
        %     for i = 1:nneurites
        %         mf_neurites(i,:) = mf_neurites_data(i).signal/femesh_neurites_mf{i}.total_volume;
        %         btpde_neurites(i,:) = btpde_neurites_data(i).signal/femesh_neurites_btpde{i}.total_volume;
                
        %     end
        % else
        %     mf_neurites = mf_neurites_data.signal/femesh_neurites_mf{1}.total_volume;
        %     btpde_neurites = btpde_neurites_data.signal/femesh_neurites_btpde{1}.total_volume;
            
        % end

        neurites_rel_error=  abs((btpde_neurites- mf_neurites)./btpde_neurites)
        neurites_abs_error=  abs(btpde_neurites- mf_neurites);

        % processes_rel_errors{i} = neurites_rel_error;
    
        namplitude = length(bvals);
        rel_error = abs((soma_btpde- soma_mf)./soma_btpde);
        abs_error  = abs(soma_btpde- soma_mf);
    
        % figure
        % plot(bvals,rel_error); hold on;
        % plot(bvals,neurites_rel_error);
        % legend(["Soma","Neurite " + string([1:nneurites])])
        % title(sprintf("Relative errors %s",cellname),'Interpreter','none');
        % grid on;
        % xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');
        % figure
        % plot(bvals,abs_error); hold on;
        % plot(bvals,neurites_abs_error);
        % legend(["Soma","Neurite " + string([1:nneurites])])
        % title(sprintf("Absolute errors %s",cellname),'Interpreter','none');
        % grid on;
        % xlabel("$b$ $\mathrm{s}/\mathrm{mm}^2$",'Interpreter','latex');
        save(sprintf('neuron_meshing_paper/microglia_output/%s_convergence_seg.mat',cellname),'rel_error','abs_error','neurites_rel_error','neurites_abs_error','bvals');   

    end
end
    

