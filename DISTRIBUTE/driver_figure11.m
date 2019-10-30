clear variables;
addpath SRC msh_files
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

%%%% Spindoctor simulation parameters.
imesh = 1; meshname_btpde{imesh} = '03b_spindle4aACC';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1; simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_dendrites_1';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_dendrites_2';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_soma';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;

imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_dendrites_1';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_dendrites_2';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_soma';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;

imesh = imesh+1; meshname_btpde{imesh} = '04b_spindle3aFI_dendrites_1';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle7aACC_dendrites_1';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;

sdeltavec = [2500,10000,10000];
ngdirvec = [90,90,90];
bdeltavec = [5000,43000,433000];
bvaluevec = [1000,4000];

nexperi = length(bdeltavec);
nb = length(bvaluevec);

%%%% Load Spindoctor simulation results
nmesh = length(meshname_btpde);
computational_time = zeros(nmesh, nb*nexperi);
number_of_FEnodes = zeros(nmesh, 1);
for imesh = 1:nmesh
    Htetgen = simul_h{imesh};
    btpde_rtol = simul_tol{imesh}(1);
    btpde_atol = simul_tol{imesh}(2);
    for ib = 1:nb
        for iexperi = 1:nexperi
            SDELTA = sdeltavec(iexperi);
            BDELTA = bdeltavec(iexperi);
            bvalue = bvaluevec(ib);
            bvecgdir_name = ['_b',num2str(bvalue),'_gdir',num2str(ngdirvec(iexperi))];
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
            
            tf = isfile([save_dir_path_spindoctor,'/',fname]);
            
            if (tf)
                call_cmd = ['load ',save_dir_path_spindoctor,'/',fname];
                disp(call_cmd);
                eval(call_cmd);
                if number_of_FEnodes(imesh)==0
                    number_of_FEnodes(imesh) = mymesh.Nnode(1);
                end
                computational_time(imesh, iexperi+nexperi*(ib-1)) = mean(ctime);
            else
                error(['cannot load ',save_dir_path_spindoctor,'/',fname]);
            end
        end
    end
end

%%%% Plot figures
marker{1} = 'r--o'; figlegend{1}='$b=1000, \delta=2.5, \Delta=5$';
marker{2} = 'r--d'; figlegend{2}='$b=1000, \delta=10, \Delta=43$';
marker{3} = 'r--*'; figlegend{3}='$b=1000, \delta=10, \Delta=433$';
marker{4} = 'k-.o'; figlegend{4}='$b=4000, \delta=2.5, \Delta=5$';
marker{5} = 'k-.d'; figlegend{5}='$b=4000, \delta=10, \Delta=43$';
marker{6} = 'k-.*'; figlegend{6}='$b=4000, \delta=10, \Delta=433$';

figure;
set(groot, 'defaultLegendInterpreter','latex');
fontsize = 14;
sorted_result = sortrows([number_of_FEnodes, computational_time], 1);
for ii=1:nb*nexperi
    plot(sorted_result(:,1)/1000, sorted_result(:, ii+1), ...
        marker{ii}, 'markersize', 10, 'LineWidth',2, 'DisplayName', figlegend{ii})
    hold on;
end
xlabel('Number of FE nodes (x1000)', 'fontsize', fontsize)
ylabel('Computational time (s)','fontsize',fontsize)
legend('location', 'northwest')
grid on;
axis tight;
set(gca,'linewidth',0.5,'fontsize',fontsize)
