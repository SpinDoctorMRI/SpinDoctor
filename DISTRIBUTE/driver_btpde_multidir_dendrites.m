clear all;
close all;
restoredefaultpath;

format short
addpath params_files/params_btpde_multidir_dendrites
addpath SRC msh_files 
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

fname_params_cells = 'params_cells_neuron.in';
fname_params_simul_domain  = 'params_simul_domain_neuron.in';
fname_params_simul_experi = 'params_simul_experi_neuron.in';

DO_PLOTS = 1;

% create cells if option 1, 2 or 3
% 1 = ellipsoids, 2 = cylinders, 3 = neuron from msh file
[params_cells,fname_cells] = create_geom(fname_params_cells);

[params_domain_geom,params_domain_pde,params_domain_femesh] ...
    = read_params_simul_domain(fname_params_simul_domain);

[DIFF_cmpts,kappa_bdys,IC_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,Ncmpt,Nboundary] ...
    = PREPARE_PDE(params_cells.ncell,params_cells.cell_shape,params_domain_geom,params_domain_pde);

fname_tetgen = [fname_cells];

%%%%%%%%%%%%% making a directory in msh
save_meshdir_path = [fname_cells,'_dir'];
tf = isfolder(save_meshdir_path);
if (~tf)
    call_cmd = ['mkdir ',save_meshdir_path];
    disp(call_cmd);
    eval([call_cmd]);
end

ns = regexp(fname_cells,'/');
save_mesh_name = [fname_cells(ns(end)+1:end),'Htetgen',num2str(params_domain_femesh.Htetgen),'msh'];
fname_tetgen = [save_meshdir_path,'/',save_mesh_name];
fname_tmp = [fname_tetgen,'.1'];

tf = isfile([fname_tmp,'.node']);
if (tf)
    fname_tetgen_femesh = fname_tmp;
else
    [fname_tetgen_femesh] = ...
    create_femesh_fromcells(params_cells,fname_cells,params_domain_geom,params_domain_femesh,fname_tetgen);
end

disp(['Reading FE mesh from:', fname_tetgen_femesh]);

if (params_cells.cell_shape == 3)
    params_cells.para_deform = [0,0]';
    [mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);
else
    [mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);
end

if (DO_PLOTS)
    PLOT_FEMESH(mymesh,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);
end

[experi_common,experi_hadc,experi_btpde] ...
    = read_params_simul_experi(fname_params_simul_experi);

[VOL_cmpts,SA_cmpts,SAu_cmpts,VOL_allcmpts,VF_cmpts,SoV_cmpts] ...
    = GET_VOL_SA(mymesh,experi_common.gdir);

ns = regexp(fname_tetgen_femesh,'/');

save_dir_name = fname_tetgen_femesh(ns(end)+1:end);

save_dir_path_spindoctor = ['saved_simul','/spindoctor/',save_dir_name];
tf = isfolder(save_dir_path_spindoctor);
if (~tf)
    call_cmd = ['mkdir ',save_dir_path_spindoctor];
    disp(call_cmd);
    eval([call_cmd]);
end

if (experi_common.ngdir_total == 1)    
    if (~isempty(experi_btpde))
        % BTPDE
        fname = ['BTPDE', '_onegdir', '_atol', num2str(experi_btpde.atol),'rtol',num2str(experi_btpde.rtol), ...
            'Htetgen',num2str(params_domain_femesh.Htetgen), '.mat'];
        if(isfile([save_dir_path_spindoctor,'/',fname]))
            call_cmd = ['load ',save_dir_path_spindoctor,'/',fname];
            eval([call_cmd]);
            disp(call_cmd);
        else
            [SOL,SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,difftime,ctime_btpde] ...
                = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,OUTPUT_MAGNETIZATION);
            [ADC_BTPDE_cmpts,ADC_BTPDE_allcmpts,ADC_allcmpts_S0] ...
                = FIT_SIGNAL(SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,experi_btpde.bvalues);
            [Sig_free,ADC_free_allcmpts] = ADCFREE(experi_btpde.bvalues,DIFF_cmpts,VOL_cmpts,IC_cmpts);
            
            call_cmd = ['save ',save_dir_path_spindoctor,'/',fname, ' ', 'experi_btpde ', 'SIG_BTPDE_allcmpts ', ...
                'Sig_free ', 'ADC_allcmpts_S0 ', 'ADC_BTPDE_allcmpts'];
            disp(['save ',save_dir_path_spindoctor,'/',fname]);
            eval([call_cmd]);
        end
        if (DO_PLOTS)
            PLOT_SIGNAL(experi_btpde.bvalues,SIG_BTPDE_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_BTPDE_allcmpts,'BTPDE');
        end
    end
    %HADC
    if (~isempty(experi_hadc))
        fname = ['HADC', '_onegdir', '_atol', num2str(experi_hadc.atol),'rtol',num2str(experi_hadc.rtol), ...
            'Htetgen',num2str(params_domain_femesh.Htetgen), '.mat'];
        if(isfile([save_dir_path_spindoctor,'/',fname]))
            call_cmd = ['load ',save_dir_path_spindoctor,'/',fname];
            eval([call_cmd]);
            disp(call_cmd);
        else
            [ADC_HADC_cmpts,ADC_HADC_allcmpts,HADC_elapsed_time] ...
            = HADC(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
            [ADC_STA_cmpts,ADC_STA_allcmpts] = STA(experi_common,DIFF_cmpts,VOL_cmpts,SAu_cmpts,IC_cmpts);
            
            call_cmd = ['save ',save_dir_path_spindoctor,'/',fname, ' ', 'ADC_HADC_cmpts ', ...
                'ADC_HADC_allcmpts ', 'DIFF_cmpts ', 'HADC_elapsed_time ', 'ADC_STA_cmpts ', 'ADC_STA_allcmpts'];
            disp(['save ',save_dir_path_spindoctor,'/',fname]);
            eval([call_cmd]);
        end
        if (DO_PLOTS)
            PLOT_ADC(ADC_HADC_cmpts,ADC_HADC_allcmpts,DIFF_cmpts,'HADC');
            PLOT_ADC(ADC_STA_cmpts,ADC_STA_allcmpts,DIFF_cmpts,'STA');
        end
    end
else
    if (~isempty(experi_btpde))
        % BTPDE
        ngdir_total = experi_btpde.ngdir_total;
        [points_gdir,graddir_index,negii] = HARDI_PTS(ngdir_total);
        
        nb = size(experi_btpde.bvalues,2);
        nexperi = size(experi_btpde.bvalues,1);
        
        for iexperi = 1:nexperi
            for ib = 1:nb
                experi_tmp = experi_btpde;
                experi_tmp.sdeltavec = experi_btpde.sdeltavec(iexperi);
                experi_tmp.bdeltavec = experi_btpde.bdeltavec(iexperi);
                experi_tmp.seqvec = experi_btpde.seqvec(iexperi);
                experi_tmp.npervec = experi_btpde.npervec(iexperi);
                experi_tmp.qvalues = experi_btpde.qvalues(iexperi,ib);
                experi_tmp.bvalues = experi_btpde.bvalues(iexperi,ib);
                
                SDELTA = experi_btpde.sdeltavec(iexperi);
                BDELTA = experi_btpde.bdeltavec(iexperi);
                bvalue = experi_btpde.bvalues(iexperi,ib);
                qvalue = experi_btpde.qvalues(iexperi,ib);
                
                var_name = ['BTPDE','_atol',num2str(experi_btpde.atol),'rtol',num2str(experi_btpde.rtol),...
                    'Htetgen',num2str(params_domain_femesh.Htetgen)];                
                [fname] = savefile_name(SDELTA,BDELTA,bvalue,ngdir_total,0);
                fname = [fname,'_',var_name,'.mat'];
                
                tf = isfile([save_dir_path_spindoctor,'/',fname]);
                
                if (tf)
                    call_cmd = ['load ',save_dir_path_spindoctor,'/',fname,...
                        ' ctime sig_cmpts sig_allcmpts'];
                    disp(call_cmd);
                    eval([call_cmd]);
                    ctime_btpde_hardi_one = ctime;
                    SIG_BTPDE_cmpts_hardi_one = sig_cmpts;
                    SIG_BTPDE_allcmpts_hardi_one = sig_allcmpts;
                else
                    [SIG_BTPDE_cmpts_hardi_one,SIG_BTPDE_allcmpts_hardi_one,ctime_btpde_hardi_one] ...
                        = SIG_BTPDE_HARDI(experi_tmp,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,...
                        points_gdir,graddir_index,negii); 
                    BTPDE.points_gdir = points_gdir;
                    BTPDE.ctime = ctime_btpde_hardi_one;
                    BTPDE.sig_cmpts = SIG_BTPDE_cmpts_hardi_one;
                    BTPDE.sig_allcmpts = SIG_BTPDE_allcmpts_hardi_one;
                    BTPDE.mymesh = mymesh;
                    BTPDE.experi_btpde = experi_tmp;
                    BTPDE.DIFF_cmpts = DIFF_cmpts;
                    BTPDE.kappa_bdys = kappa_bdys;
                    BTPDE.params_cells = params_cells;
                    BTPDE.params_domain_geom = params_domain_geom;
                    BTPDE.params_domain_femesh = params_domain_femesh;
                    BTPDE.params_domain_pde = params_domain_pde;
                    BTPDE.VOL_cmpts = VOL_cmpts;
                    BTPDE.SA_cmpts = SA_cmpts;
                    BTPDE.IC_cmpts = IC_cmpts;
                    BTPDE.IN_cmpts_index = IN_cmpts_index;
                    BTPDE.OUT_cmpts_index = OUT_cmpts_index;
                    BTPDE.ECS_cmpts_index = ECS_cmpts_index;
                    BTPDE.cmpts_bdys_mat= cmpts_bdys_mat;
                    call_cmd = ['save ',save_dir_path_spindoctor,'/',fname, ' ', '-struct', ' ', 'BTPDE'];
                    disp(call_cmd);
                    eval([call_cmd]);
                end
                
                SIG_BTPDE_cmpts_hardi(:,:,iexperi,ib) = SIG_BTPDE_cmpts_hardi_one;
                SIG_BTPDE_allcmpts_hardi(:,iexperi,ib) = SIG_BTPDE_allcmpts_hardi_one;
                ctime_btpde_hardi(:,iexperi,ib) = ctime_btpde_hardi_one;
            end
        end
        
        if (DO_PLOTS)
            cmpts_use = unique([1,Ncmpt]);
            experi_use = unique([1,nexperi]);
            bvec_use = unique([2,nb]);

            for iexperi = experi_use
                for ib = bvec_use
                    S0 = sum((IC_cmpts.*VOL_cmpts));
                    bv = experi_btpde.bvalues(iexperi,ib);
                    title_str = ['BTPDE. All Cmpts. ','Experi ',...
                        num2str(iexperi),', b= ',num2str(bv)];
                    PLOT_HARDI_PT(points_gdir,...
                        squeeze(real(SIG_BTPDE_allcmpts_hardi(:,iexperi,ib)))/S0,title_str);
                end
            end
        end
    end
    
    %HADC
    if (~isempty(experi_hadc))
        fname = ['HADC_HARDI','_gdir', num2str(experi_hadc.ngdir_total), '_atol',num2str(experi_hadc.atol),'rtol',num2str(experi_hadc.rtol),...
            'Htetgen',num2str(params_domain_femesh.Htetgen), '.mat'];
        if (isfile([save_dir_path_spindoctor,'/',fname]))
            call_cmd = ['load ',save_dir_path_spindoctor,'/',fname];
            eval([call_cmd]);
            disp(call_cmd);
        else
            ngdir_total = experi_hadc.ngdir_total;
            [points_gdir,graddir_index,negii] = HARDI_PTS(ngdir_total);
            
            [ADC_HADC_cmpts_hardi,ADC_HADC_allcmpts_hardi,ctime_hadc_hardi] ...
                = HADC_HARDI(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts,...
                points_gdir,graddir_index,negii);
        
            HARDI.experi_hadc = experi_hadc;
            HARDI.points_gdir = points_gdir;
            HARDI.ADC_HADC_cmpts_alldir = ADC_HADC_cmpts_hardi;
            HARDI.ADC_HADC_allcmpts_alldir = ADC_HADC_allcmpts_hardi;
            call_cmd = ['save ',save_dir_path_spindoctor,'/',fname, ' ', '-struct', ' ', 'HARDI'];
            disp(call_cmd);
            eval([call_cmd]);
        end
        
        if(DO_PLOTS)
            nexperi = length(experi_hadc.sdeltavec);
            for iexperi = 1:nexperi
                PLOT_HARDI_PT(HARDI.points_gdir,HARDI.ADC_HADC_allcmpts_alldir(:,iexperi),['HADC\_HARDI experiment', num2str(iexperi)]);
            end
        end
    end
end
