function driver_get_morez_errors(meshname, lsval,hval ,segment_cell,direc)    

nh = length(hval);
nls = length(lsval);
segment_cell = str2num(segment_cell);

if ~isdir(direc)
    mkdir(direc)
end

[mesh_path,cellname,~] = fileparts(meshname);
swc_file=sprintf("swc_files/%s.swc",cellname);
addpath(genpath('src'));
addpath(genpath('setups'));
setup_morez_ref_sol;

setup_r = setup;
setup_r.name = string(meshname);

setup_r.geometry.tetgen_options = sprintf("-pq1.2a%.1fVCn",0.1);
[setup_r,femesh_btpde,~,~] = prepare_simulation(setup_r);
savepath_r=fullfile('saved_simul',sprintf('%s_tet%s',cellname,setup_r.geometry.tetgen_options),'cell');
btpde_cell = load_btpde(setup_r,savepath_r,false);

if segment_cell ==1
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

rel_error = cell(nh,nls);
if segment_cell ==1
    rel_error_soma =  cell(nh,nls);
    rel_error_dendrites = cell(ndendrites,nh,nlss);
end
styles = ["-";"--";":"];
colors = ['g';'b';'r'];

for i = 1:nh
    for j = 1:nls
        clear setup;
        setup_morez_ref_sol;
        h = hval(i); ls = lsval(j);
        setup.name = string(meshname);
        setup.geometry.tetgen_options = sprintf("-pq1.2a%.1fVCn",h);
        setup.mf.length_scale = ls;
        [setup,femesh,~,~] = prepare_simulation(setup);
        savepath=fullfile('saved_simul',sprintf('%s_tet%s',cellname,setup.geometry.tetgen_options),'cell');
        mf_cell = load_mf(setup,savepath,false);
        rel_error = abs(real(mf_cell.signal - btpde_cell.signal))./abs(real(btpde_cell.signal))

        fig_cell = figure(1); hold on;
        name = sprintf("h= %.1f, ls =%.4f",h,ls);
        plot(rel_error,'Color',colors(j),'LineStyle',styles(i),'DisplayName',name);

        if segment_cell == 1
            tetgen_path=sprintf('%s/%s_ply_dir/%s_%s_tet%s_mesh.1',mesh_path,cellname,cellname,setup.geometry.ecs_shape,setup.geometry.tetgen_options);
            [~,femesh_dendrites] = segment_femesh(femesh,swc_file,tetgen_path);
            permutation = align_dendrites(femesh_dendrites_btpde,femesh_dendrites);
            savepath_soma=fullfile('saved_simul',sprintf('%s_tet%s',cellname,setup.geometry.tetgen_options),'soma');
            mf_soma = load_mf(setup,savepath_soma,false);
            rel_error_dendrites = cell(ndendrites,1);
            for k =1:ndendrites
                savepath_dendrite =fullfile('saved_simul',sprintf('%s_tet%s/dendrite_%d',cellname,setup.geometry.tetgen_options,permutation(k)));
                mf_dendrite = load_mf(setup,savepath_dendrite,false);
                rel_error_dendrites{k} = abs(real(mf_dendrite.signal - btpde_dendrite{k}.signal))./abs(real(btpde_dendrite{k}.signal))
            end
            rel_error_soma = abs(real(mf_soma.signal - btpde_soma.signal))./abs(real(btpde_soma.signal))
            
            fig_soma = figure(2); hold on;
            plot(rel_error_soma,'Color',colors(j),'LineStyle',styles(i),'DisplayName',name);
            
            for k = 1:ndendrites
                fig_dend = figure(2+k);hold on;
                plot(rel_error_dendrites{k},'Color',colors(j),'LineStyle',styles(i),'DisplayName',name);
            end
        end
        
        


    end
end
fig_cell = figure(1);grid on;xlabel('Sequence');ylabel('Relative error');legend('Location','eastoutside');
saveas(fig_cell,sprintf('%s/%s_test_morez_cell.fig',direc,cellname))

if segment_cell==1
    fig_soma = figure(2);grid on;xlabel('Sequence');ylabel('Relative error');legend('Location','eastoutside');
    saveas(fig_soma,sprintf('%s/%s_test_morez_soma.fig',direc,cellname))

    for k = 1:ndendrites
        fig_dend = figure(2+k);grid on;xlabel('Sequence');ylabel('Relative error');legend('Location','eastoutside');
        saveas(fig_dend,sprintf('%s/%s_test_morez_dendrite_%d.fig',direc,cellname,k))
    end
end