clear all;
addpath SRC msh_files
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

for jj=1:2
    if jj==1
        % Figure 10A
        imesh = 1; meshname_btpde{imesh} = '03a_spindle2aFI';
        simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = -1; simul_2d{imesh}=1;
        imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_dendrites_1';
        simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = -1;simul_2d{imesh}=1;
        imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_dendrites_2';
        simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = -1;simul_2d{imesh}=1;
        sdeltavec = 10000;
        bdeltavec = 43000;
        bvaluevec = [60000, 40000, 20000, 12000, 10000, 8000, 7000, 6000, 4000, 2500];
        S_avg_dendrite1 = @(b) -0.0221+24.12.*b.^(-1/2); S_avg_dendrite1_string = '$S=-0.0221+24.12\cdot b^{-1/2}$';
        S_avg_dendrite2 = @(b) -0.0154+24.96.*b.^(-1/2); S_avg_dendrite2_string = '$S=-0.0154+24.96\cdot b^{-1/2}$';    
    else
        % Figure 10B
        imesh = 1; meshname_btpde{imesh} = '03b_spindle4aACC';
        simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = -1; simul_2d{imesh}=1;
        imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_dendrites_1';
        simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = -1;simul_2d{imesh}=1;
        imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_dendrites_2';
        simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = -1;simul_2d{imesh}=1;
        sdeltavec = 10000;
        bdeltavec = 43000;
        bvaluevec = [60000, 40000, 20000, 12000, 10000, 8000, 7000, 6000, 4000, 2500];
        S_avg_dendrite1 = @(b) 0.0028+24.05.*b.^(-1/2); S_avg_dendrite1_string = '$S=0.0028+24.05 \cdot b^{-1/2}$';
        S_avg_dendrite2 = @(b) 0.0322+23.94.*b.^(-1/2); S_avg_dendrite2_string = '$S=0.0322+23.94 \cdot b^{-1/2}$';
    end
    
    nexperi = length(bdeltavec);
    nb = length(bvaluevec);
    nmesh = length(meshname_btpde);
    
    sig_avg_b = zeros(nb+1, nmesh);
    for imesh = 1:nmesh
        Htetgen = simul_h{imesh};
        btpde_rtol = simul_tol{imesh}(1);
        btpde_atol = simul_tol{imesh}(2);
        for iexperi = 1:nexperi
            for ib = 1:nb+1
                SDELTA = sdeltavec(iexperi);
                BDELTA = bdeltavec(iexperi);
                if ib < nb+1
                    bvalue = bvaluevec(ib);
                else
                    bvalue = 0;
                end
                bvecgdir_name = ['_b',num2str(bvalue),'_gdir30'];
                experi_name = ['d',num2str(SDELTA),'_D',num2str(BDELTA),bvecgdir_name];
                save_dir_name_spindoctor = [meshname_btpde{imesh},'Htetgen',num2str(Htetgen),'msh.1'];
                save_dir_path_spindoctor = ['saved_simul','/spindoctor/',save_dir_name_spindoctor];
                var_name = ['BTPDE','_atol',num2str(btpde_atol),'rtol',num2str(btpde_rtol),...
                    'Htetgen',num2str(Htetgen), '.mat'];
                if (simul_2d{imesh} == 0)
                    fname = [experi_name,'_cmpt0','_',var_name];
                else
                    fname = [experi_name,'_cmpt0','gdir2D_',var_name];
                end
                
                if (isfile([save_dir_path_spindoctor,'/',fname]))
                    call_cmd = ['load ',save_dir_path_spindoctor,'/',fname];
                    disp(call_cmd);
                    eval(call_cmd);
                    sig_avg_b(ib, imesh) = mean(real(sig_allcmpts));
                else
                    error(['cannot load ',save_dir_path_spindoctor,'/',fname]);
                end
            end
        end
    end
    
    sig_avg_normalized = sig_avg_b(1:end-1 ,:) ./ sig_avg_b(end, :);
    inverse_sqrtb = 1./sqrt(bvaluevec);
    
    figure;
    figlegend{1}='Neuron'; marker{1} = 'co';
    figlegend{2}='Dendrite 1'; marker{2} = 'rd';
    figlegend{3}='Dendrite 2'; marker{3} = 'k>';
    set(groot, 'defaultLegendInterpreter','latex');
    fontsize = 14;
    for ii=1:nmesh
        plot(inverse_sqrtb, sig_avg_normalized(:, ii), ...
            marker{ii}, 'markersize', 10, 'LineWidth',2, 'DisplayName', figlegend{ii})
        hold on;
    end
    plot(1./sqrt(linspace(2500,60000,1000)), S_avg_dendrite1(linspace(2500,60000,1000)), ...
        'b-', 'LineWidth', 2, 'DisplayName', S_avg_dendrite1_string)
    plot(1./sqrt(linspace(2500,60000,1000)), S_avg_dendrite2(linspace(2500,60000,1000)), ...
        'b--', 'LineWidth', 2, 'DisplayName', S_avg_dendrite2_string)
    xlabel('$1/\sqrt{b}(mm \cdot s^{-1/2})$', 'fontsize', fontsize, 'Interpreter','latex')
    ylabel('$S_{ave}(b)$','fontsize',fontsize, 'Interpreter','latex')
    title(meshname_btpde{1},'fontsize',fontsize,'Interpreter', 'none')
    legend('location', 'northwest')
    grid on;
    axis tight;
    set(gca,'linewidth',0.5,'fontsize',fontsize)
end
