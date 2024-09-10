function results = driver_get_PGSE_errors(meshname,only_cell,setup_file,direc,lsval,hval) 
% DRIVER_GET_PGSE_ERRORS Records 

hval = str2num(hval);
lsval = str2num(lsval);

set(groot,'defaultLineLineWidth',3.0)
nh = length(hval);
nls = length(lsval);
only_cell = str2num(only_cell);

if ~isdir(direc)
    mkdir(direc)
end

[mesh_path,cellname,~] = fileparts(meshname);
swc_file=sprintf("swc_files/%s.swc",cellname);
addpath(genpath('src'));
addpath(genpath('setups'));
fprintf("Running %s.m\n",setup_file)
run(sprintf("%s.m",setup_file));
setup_r = setup;
setup_r.name = string(meshname);

setup_r.geometry.tetgen_options = "-pq1.2a0.1VCn";
[setup_r,femesh_btpde,~,~] = prepare_simulation(setup_r);
savepath_r=fullfile('saved_simul',sprintf('%s_tet%s',cellname,setup_r.geometry.tetgen_options),'cell');
btpde_cell = load_btpde(setup_r,savepath_r,false);

if only_cell ~=1
    savepath_soma_r=fullfile('saved_simul',sprintf('%s_tet%s',cellname,setup_r.geometry.tetgen_options),'soma');
    btpde_soma = load_btpde(setup_r,savepath_soma_r,false);

    tetgen_path_btpde=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup_r.geometry.ecs_shape,setup_r.geometry.tetgen_options);
    [~,femesh_dendrites_btpde] = segment_femesh(femesh_btpde,swc_file,tetgen_path_btpde);
    ndendrites = length(femesh_dendrites_btpde);
    for k =1:ndendrites
        savepath_dendrite =fullfile('saved_simul',sprintf('%s_tet%s/dendrite_%d',cellname,setup_r.geometry.tetgen_options,k));
        btpde_dendrite{k} = load_btpde(setup_r,savepath_dendrite,false);
    end
end

% rel_error = cell(nh,nls);
% if only_cell ~=1
%     rel_error_soma =  cell(nh,nls);
%     rel_error_dendrites = cell(ndendrites,nh,nls);
% end
styles = ["-";"--";":"];
colors = ['k';'r'];
if only_cell ~= 1
    cm = colormap(parula(ndendrites));  
end
max_errors = zeros(nh,nls);
results = cell(nh,nls);
for i = 1:nh
    for j = 1:nls
        try
        clear setup;
        run(sprintf("%s.m",setup_file));
        output = struct;
        h = hval(i); ls = lsval(j);
        output.h = h; output.ls = ls;
        setup.name = string(meshname);
        setup.geometry.tetgen_options = sprintf("-pq1.2a%.1fVCn",h);%sprintf("-pq1.2a%.1fVCn",h);
        setup.mf.length_scale = ls;
        [setup,femesh,~,~] = prepare_simulation(setup);
        savepath=fullfile('saved_simul',sprintf('%s_tet%s',cellname,setup.geometry.tetgen_options),'cell');
        mf_cell = load_mf(setup,savepath,false);
        rel_error = abs(real(mf_cell.signal - btpde_cell.signal))./abs(real(btpde_cell.signal));
        output.cell = rel_error;
        fig_cell = figure(i); hold on;
        name = sprintf("h= %.1f, ls =%.4f, cell",h,ls);
        % plot(setup.gradient.bvalues,rel_error,'Color',colors(i),'LineStyle',styles(j),'DisplayName',name);
        plot(setup.gradient.bvalues,rel_error,'LineStyle',styles(j),'DisplayName',name,'Color',colors(1),'Marker','x');

        if only_cell ~= 1
            tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
            [~,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);
            permutation = align_dendrites(femesh_dendrites_btpde,femesh_dendrites);
            savepath_soma=fullfile('saved_simul',sprintf('%s_tet%s',cellname,setup.geometry.tetgen_options),'soma');
            mf_soma = load_mf(setup,savepath_soma,false);
            rel_error_dendrites = cell(ndendrites,1);
            for k =1:ndendrites
                savepath_dendrite =fullfile('saved_simul',sprintf('%s_tet%s/dendrite_%d',cellname,setup.geometry.tetgen_options,permutation(k)));
                mf_dendrite = load_mf(setup,savepath_dendrite,false);
                rel_error_dendrites{k} = abs(real(mf_dendrite.signal - btpde_dendrite{k}.signal))./abs(real(btpde_dendrite{k}.signal));
            end
            rel_error_soma = abs(real(mf_soma.signal - btpde_soma.signal))./abs(real(btpde_soma.signal));
            output.soma = rel_error_soma;
            output.dendrites = rel_error_dendrites;
            hold on;
            % fig_soma = figure(2); hold on;
            name = sprintf("h= %.1f, ls =%.4f, soma",h,ls);
                
            % plot(setup.gradient.bvalues,rel_error_soma,'Color',colors(i),'LineStyle',styles(j),'DisplayName',name);
            plot(setup.gradient.bvalues,rel_error_soma,'LineStyle',styles(j),'DisplayName',name,'Color',colors(2),'Marker','x');

            for k = 1:ndendrites
                % fig_dend = figure(2+k);hold on;
                name = sprintf("h= %.1f, ls =%.4f, dendrite_%d",h,ls,k);
                % plot(setup.gradient.bvalues,rel_error_dendrites{k},'Color',colors(i),'LineStyle',styles(j),'DisplayName',name);
                plot(setup.gradient.bvalues,rel_error_dendrites{k},'LineStyle',styles(j),'DisplayName',name,'Color',cm(k,:),'Marker','x');
            end
            max_errors(i,j) = max([[rel_error_dendrites{:}],rel_error,rel_error_soma]);

        else
            max_errors(i,j) = max(rel_error);
        end
        results{i,j} = output;
        catch
            warning('Data not loading for h=%f,ls=%f',h,ls);
            max_errors(i,j) = NaN;
        end

    end
    end
for i =1:nh
fig_cell = figure(i);grid on;xlabel('b-values');ylabel('Relative error');legend('Location','eastoutside','Interpreter','None');
set(gca, 'YScale', 'log');
saveas(fig_cell,sprintf('%s/%s_test_pgse_h%.1f.png',direc,cellname,hval(i)))
end
for i =1:nh
    for j = 1:nls
        fprintf("Max relative error for h = %.1f, ls = %.1f is %f\n",hval(i),lsval(j),max_errors(i,j));
    end
end
