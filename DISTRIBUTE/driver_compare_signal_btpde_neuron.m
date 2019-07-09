clear all;
addpath SRC msh_files
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

%%%% choose which experiments to display
%%%% Spindoctor simulation parameters.

imesh = 0;

% Figure 1
% imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC';
% simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1; simul_2d{imesh}=1;
% imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_dendrites_1';
% simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
% imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_dendrites_2';
% simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
% imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle4aACC_soma';
% %simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = 0.1;simul_2d{imesh}=1;
% simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
% sdeltavec = [2500,10000,10000];
% ngdirvec = [90,90,90];
% bdeltavec = [5000,43000,433000];
% bvaluevec = [1000,4000];

% Figure 2
imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_dendrites_1';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_dendrites_2';
simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
imesh = imesh+1; meshname_btpde{imesh} = '03a_spindle2aFI_soma';
simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = 0.1;simul_2d{imesh}=1;
sdeltavec = [2500,10000,10000];
ngdirvec = [90,90,90];
bdeltavec = [5000,43000,433000];
bvaluevec = [1000,4000];

% Figure 3
% imesh = imesh+1; meshname_btpde{imesh} = '04b_spindle3aFI_dendrites_1';
% simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
% imesh = imesh+1; meshname_btpde{imesh} = '03b_spindle7aACC_dendrites_1';
% simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=1;
% sdeltavec = [2500,10000,10000];
% ngdirvec = [90,90,90];
% bdeltavec = [5000,43000,433000];
% bvaluevec = [1000,4000];

% Figure 4
% imesh = imesh+1; meshname_btpde{imesh} = '25o_pyramidal18aFI';
% simul_tol{imesh} = [1e-2,1e-4]; simul_h{imesh} = -1;simul_2d{imesh}=0;
% imesh = imesh+1; meshname_btpde{imesh} = '25o_pyramidal18aFI';
% simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = -1;simul_2d{imesh}=0;
% imesh = imesh+1; meshname_btpde{imesh} = '25o_pyramidal18aFI';
% simul_tol{imesh} = [1e-3,1e-6]; simul_h{imesh} = 0.2;simul_2d{imesh}=0;
% sdeltavec = [10000];
% ngdirvec = [10];
% bdeltavec = [13000];
% bvaluevec = [1000,4000];

nexperi = length(bdeltavec);
nb = length(bvaluevec);

%%%%%%%%%%%%%%%%
%%%% Load Spindoctor simulations.
nmesh = length(meshname_btpde);

for imesh = 1:nmesh
    
    Htetgen = simul_h{imesh};
    btpde_rtol = simul_tol{imesh}(1);
    btpde_atol = simul_tol{imesh}(2);
    for iexperi = 1:nexperi
        for ib = 1:nb
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
                bt_sig_allcmpts{imesh}{iexperi}{ib} = sig_allcmpts;
                bt_ctime{imesh}{iexperi}{ib} = ctime;
                bt_bvec{imesh}{iexperi}{ib} = experi_btpde.bvalues;
                bt_rtol{imesh} = experi_btpde.rtol;
                bt_atol{imesh} = experi_btpde.atol;
                bt_htetgen{imesh} = params_domain_femesh.Htetgen;
                bt_ngdir{imesh} = experi_btpde.ngdir_total;
                bt_vol{imesh} = sum(VOL_cmpts);
                if (~exist('points_gdir'))
                    ngdir_total = experi_btpde.ngdir_total;
                    [points_gdir,graddir_index,negii] = HARDI_PTS(ngdir_total);
                end
                bt_points_gdir{imesh}{iexperi}{ib} = points_gdir;
                
                %     prt_str{imesh} = ['rtol=',num2str(bt_rtol{imesh}),'Htetgen',num2str(bt_htetgen{imesh})];
                prt_str{imesh}{iexperi}{ib} ...
                    = [meshname_btpde{imesh},',\delta=',num2str(SDELTA),',\Delta=',num2str(BDELTA),...
                    ',b=',num2str(bt_bvec{imesh}{iexperi}{ib})];
                mymesh_save{imesh} = mymesh;
            else
                disp(['cannot load ',save_dir_path_spindoctor,'/',fname]);
                bt_sig_allcmpts{imesh}{iexperi}{ib} = [];
            end
            
        end
    end
end

% ngdir_total = experi_btpde.ngdir_total;
% [points_gdir,graddir_index,negii] = HARDI_PTS(ngdir_total);

markercell_vec{1}='d';
markercell_vec{2}='o';
markercell_vec{3}='p';
markercell_vec{4}='>';
markercell_vec{5}='<';
markercell_vec{6}='x';

colorcell_vec{1} = 'r';
colorcell_vec{2} = 'c';
colorcell_vec{3} = 'b';
colorcell_vec{4} = 'g';

linestylecell_vec{1}='--';
linestylecell_vec{2}='-';
linestylecell_vec{3}='-.';
linestylecell_vec{4}=':';

% Plot figure

if (1 == 0)
    angle_max = [];
    angle_min = [];
    for imesh = 1:nmesh
        figure; hold on
        iplot = 0;
        for ib = 1:nb
            for iexperi = 1:nexperi
                if (~isempty(bt_sig_allcmpts{imesh}{iexperi}{ib}))
                    ii = find(abs(bt_points_gdir{imesh}{iexperi}{ib}(:,3) - 0)<=1e-6);
                    iplot = iplot + 1;
                    bvec = bt_bvec{imesh}{iexperi}{ib};
                    xvec = atan2(bt_points_gdir{imesh}{iexperi}{ib}(:,2),bt_points_gdir{imesh}{iexperi}{ib}(:,1));
                    angle_all{imesh}{iexperi}{ib} = xvec(ii);
                    xvec = atan(bt_points_gdir{imesh}{iexperi}{ib}(:,2)./bt_points_gdir{imesh}{iexperi}{ib}(:,1));
                    yvec = real(bt_sig_allcmpts{imesh}{iexperi}{ib}/bt_vol{imesh});
                    xx = xvec(ii);
                    yy = yvec(ii);
                    
                    
                    [y,imax] = max(yy);
                    [y,imin] = min(yy);
                    xx = atan(sin(xx)./cos(xx));
                    angle_max(imesh,iexperi,ib) = xx(imax);
                    angle_min(imesh,iexperi,ib) = xx(imin);
                    
                    [xval,jj] = sort(xx);
                    h = plot((xx(jj)),(yy(jj)),...
                        [colorcell_vec{mod(iexperi-1,4)+1},'.']);%markercell_vec{mod(ib-1,4)+1}]);
                    set(h,'MarkerSize', 8, 'LineWidth',1);
                    legend_vec{iplot} = [prt_str{imesh}{iexperi}{ib}];
                end
            end
        end
        xlabel('gdir');
        ylabel('Signal');
        legend(legend_vec{1:iplot});
        %legend('Location','southwest');
        set(gca,'ylim',[0,1]);
        grid on;
    end
end



for imesh = 1:nmesh
    iplot = 0;
    
    xmin=min(mymesh_save{imesh}.Pts_cmpt_reorder{1}(1,:));
    xmax=max(mymesh_save{imesh}.Pts_cmpt_reorder{1}(1,:));
    ymin=min(mymesh_save{imesh}.Pts_cmpt_reorder{1}(2,:));
    ymax=max(mymesh_save{imesh}.Pts_cmpt_reorder{1}(2,:));
    %rlen = max(abs([xmax,xmin,ymax,ymin]));
    rlen = 200;
    
    mymesh_plot = mymesh_save{imesh};
    mymesh_plot.Pts_cmpt_reorder{1} = mymesh_plot.Pts_cmpt_reorder{1}/rlen*bt_vol{imesh};
    
    PLOT_FEMESH(mymesh_plot,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);
    title([meshname_btpde{imesh},',VOL=',num2str(bt_vol{imesh}),], 'Interpreter', 'none');
    view([0,0,1]); hold on;
    legend_vec{1} = 'cell';
    iplot = 1;   
    
    %figure; hold on;
    rlen = bt_vol{imesh}; %sqrt((xlen/2)^2+(ylen/2)^2);
    

    for ib = 1:nb
        for iexperi = 1:nexperi
            if (~isempty(bt_sig_allcmpts{imesh}{iexperi}{ib}))
                ii = find(abs(bt_points_gdir{imesh}{iexperi}{ib}(:,3) - 0)<=1e-6);
                xvec = atan(bt_points_gdir{imesh}{iexperi}{ib}(:,2)./bt_points_gdir{imesh}{iexperi}{ib}(:,1));
                yvec = real(bt_sig_allcmpts{imesh}{iexperi}{ib}); %/bt_vol{imesh};
                xx = xvec(ii);
                yy = yvec(ii);
                
                [xval,jj] = sort(xx);
                angle_use = xx(jj);
                
                cvx = cos(angle_use).*yy(jj);
                cvy = sin(angle_use).*yy(jj);
                xpt=[cvx;-cvx];
                ypt=[cvy;-cvy];
                h=plot(xpt,ypt,[colorcell_vec{mod(iexperi-1,4)+1},...
                    linestylecell_vec{mod(ib-1,4)+1}]);
                set(h,'markersize',8,'linewidth',2);
                iplot = iplot+1;
                legend_vec{iplot} = ['\Delta=',num2str(bdeltavec(iexperi)/1000),', b=',num2str(bvaluevec(ib))];
                
            end
        end
    end
    legend(legend_vec{1:iplot},'Location','northwest');
    set(gca,'xlim',[-rlen,rlen]);
    set(gca,'ylim',[-rlen,rlen]);
    axis equal;
    grid on;
    set(gca,'xlim',[-rlen,rlen]);
    set(gca,'ylim',[-rlen,rlen]);
end

if (1 == 1)
    % Plot figure
    figure; hold on
    iplot = 0;   
    for imesh = 1:nmesh-1       
        for iexperi = 1:nexperi
            for ib = 1:nb
                yvec_ref = real(bt_sig_allcmpts{end}{iexperi}{ib}./bt_vol{end});
                iplot = iplot + 1;
                bvec = bt_bvec{imesh}{iexperi}{ib};
               
                %xvec = atan(bt_points_gdir{imesh}{iexperi}{ib}(:,2)./bt_points_gdir{imesh}{iexperi}{ib}(:,1));
                yvec = real(bt_sig_allcmpts{imesh}{iexperi}{ib}./bt_vol{imesh});
                xvec = 1:length(yvec);
               
                h = plot((xvec),abs((yvec-yvec_ref)./yvec_ref)*100,[colorcell_vec{mod(imesh-1,4)+1},...
                    linestylecell_vec{mod(ib-1,4)+1}]);
                set(h,'markersize',8,'linewidth',2);

                
                legend_vec{iplot} = ...
                    ['Htetgen=',num2str(simul_h{imesh}), ', ODE tol=[',num2str(simul_tol{imesh}),...
                    '], \Delta=',num2str(bdeltavec(iexperi)/1000),', b=',num2str(bvaluevec(ib))];
            end
        end
    end
    
    xlabel('gdir number');
    ylabel('Relative signal difference (%)');
    legend(legend_vec{1:iplot});
    % legend('Location','southwest');
    %set(gca,'ylim',[0,1]);
    grid on;    
end

