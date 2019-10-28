%%%% The DISTRIBUTE folder contains a commented general purpose driver 
%%%% called driver_spindoctor_neuronmodule_commented.m. 
%%%% It is highly recommended to read this driver to understand the workflow 
%%%% of SpinDoctor and the Neuron Module. 
%%%% driver_spindoctor_neuronmodule_commented.m does not use saved simulation data, 
%%%% all simulations are run from scratch.
%%%% The user is advised to read the latest version 
%%%% from \url{https://github.com/jingrebeccali/SpinDoctor} 

clear all;
close all;
restoredefaultpath;

%%%% User needs to choose to see some of the typical plots or not.
DO_PLOTS = true;

%%%% User needs to choose whether to output the solution of the BTPDE, 
%%%% in the case only one diffusion direction is simulated.  
%%%% If the finite elements mesh is large, it is
%%%% recommended not to output the magnetization.  
%%%% If the flag is set to false, the function PLOT_PDESOLUTION cannot be
%%%% called.
%%%% For multiple diffusion directions (HARDI), the magnetization is never
%%%% outputted and the flag below is not taken into account.
OUTPUT_MAGNETIZATION = true;

%%%% User needs to define the directory where the 3 user provided input files are kept.
user_inputfiles_dir = 'params_files/params_spindoctor_neuronmodule_commented';

%%%% User needs to define the names of the 3 user provided input files
%%%% found in user_inputfiles_dir
fname_params_cells = 'params_cells_neuron.in';
fname_params_simul_domain  = 'params_simul_domain_neuron.in';
fname_params_simul_experi = 'params_simul_experi_neuron.in';

%%%% We add the path to the directory of the 3 user provided input files.
eval(['addpath ', user_inputfiles_dir]);

%%%% We add the path to the top level source code directory.
addpath SRC
%%%% We add the path to the deeper level source code directories.
%%%% We note the functions that need to be called in a typical simulation
%%%% situation are contained in the top level directory 'SRC'.
%%%% The less often used functions are kept in the deeper levels.
addpath SRC/PDE SRC/DMRI SRC/FEM SRC/GEOM SRC/TETGEN

%%%% We read the params_cells input file and create the geometrical configuration
%%%% In the first line of the params_cells input file, we expect a number
%%%% between 1 and 3;
%%%% 1 = ellipsoids, 2 = cylinders, 3 = neuron from msh file;
[params_cells,fname_cells] = create_geom(fname_params_cells);

%%%% We read the params_simul_domain input file;
[params_domain_geom,params_domain_pde,params_domain_femesh] ...
    = read_params_simul_domain(fname_params_simul_domain);

%%%% We set up the PDE model in the geometrical compartments.
[DIFF_cmpts,kappa_bdys,IC_cmpts,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,Ncmpt,Nboundary] ...
    = PREPARE_PDE(params_cells.ncell,params_cells.cell_shape,params_domain_geom,params_domain_pde);

%%%% We make a directory to store the finite elements mesh for this simulation;
save_meshdir_path = [fname_cells,'_dir'];
tf = isfolder(save_meshdir_path);
if (~tf)
    call_cmd = ['mkdir ',save_meshdir_path];
    disp(call_cmd);
    eval([call_cmd]);
end
%%%% We use an existing finite elements mesh or we create a new finite
%%%% elements mesh.
%%%% The name of the finite elements mesh is stored in the string fname_tetgen_femesh
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


%%%% We do not deform the finite elements mesh if in the Neuron module.
%%%% We do deform the finite elements mesh if in SpinDoctor White Matter
%%%% mode.
%%%% The finite elements mesh is stored in mymesh
if (params_cells.cell_shape == 3)
    params_cells.para_deform = [0,0]';
    [mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);
else
    [mymesh,cmpts_bdys_mat] = read_tetgen(fname_tetgen_femesh,params_cells.para_deform,Ncmpt,Nboundary);
end

%%%% We plot the finite elements mesh
if (DO_PLOTS)
    PLOT_FEMESH(mymesh,OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index);
end

%%%% We read the params_simul_experi input file;
[experi_common,experi_hadc,experi_btpde] ...
    = read_params_simul_experi(fname_params_simul_experi);

%%%% We get the volume and the surface area quantities from the mesh
[VOL_cmpts,SA_cmpts,SAu_cmpts,VOL_allcmpts,VF_cmpts,SoV_cmpts] ...
    = GET_VOL_SA(mymesh,experi_common.gdir);

nexperi = length(experi_common.sdeltavec);
 
%%%% We run the simulation for one diffusion-encoding direction
if (experi_common.ngdir_total == 1)
    if (~isempty(experi_btpde))        
        nb = size(experi_btpde.bvalues,2);         
        %%%% We solve the BTPDE        
        [SOL,SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,difftime,ctime_btpde] ...
            = BTPDE(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,OUTPUT_MAGNETIZATION);
        %%%% We fit the ADC of the signal
        [ADC_BTPDE_cmpts,ADC_BTPDE_allcmpts,ADC_allcmpts_S0] ...
            = FIT_SIGNAL(SIG_BTPDE_cmpts,SIG_BTPDE_allcmpts,experi_btpde.bvalues);
        %%%% We obtain the signal of free diffusion
        [Sig_free,ADC_free_allcmpts] = ADCFREE(experi_btpde.bvalues,DIFF_cmpts,VOL_cmpts,IC_cmpts);
        %%%% We plot the BTPDE signal, the S0*exp(-ADC*b) curve, and the
        %%%% free diffusion curves together.
        if (DO_PLOTS)
            PLOT_SIGNAL(experi_btpde.bvalues,SIG_BTPDE_allcmpts,Sig_free,ADC_allcmpts_S0,ADC_BTPDE_allcmpts,'BTPDE');
        end
        %%%% We plot the magnetization solution of the BTPDE if the user
        %%%% has set the flag OUTPUT_MAGNETIZATION to true
        if (DO_PLOTS & OUTPUT_MAGNETIZATION)
            for iexperi = 1:nexperi
                for ib = 1:nb
                    bv = experi_btpde.bvalues(iexperi,ib);
                    title_str = ['BTPDE magnetization. ','Experi ',...
                        num2str(iexperi),', b= ',num2str(bv)]; 
                    PLOT_PDESOLUTION(mymesh,SOL{iexperi}{ib},OUT_cmpts_index,ECS_cmpts_index,IN_cmpts_index,title_str);
                end
            end
        end
    end    
    if (~isempty(experi_hadc))      
        %%%% We solve the HADC model
        [ADC_HADC_cmpts,ADC_HADC_allcmpts,ctime_btpde] ...
            = HADC(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts);
        %%%% We find the short time approximation (STA) of the effective
        %%%% diffusion.
        [ADC_STA_cmpts,ADC_STA_allcmpts] = STA(experi_common,DIFF_cmpts,VOL_cmpts,SAu_cmpts,IC_cmpts);   
        %%%% We plot the HADC and the STA
        if (DO_PLOTS)
            PLOT_ADC(ADC_HADC_cmpts,ADC_HADC_allcmpts,DIFF_cmpts,'HADC');            
            PLOT_ADC(ADC_STA_cmpts,ADC_STA_allcmpts,DIFF_cmpts,'STA');
        end
    end
end
%%%% We run the simulation for many diffusion-encoding directions uniformly
%%%% distributed in the unit sphere in 3 dimensions.
if (experi_common.ngdir_total > 1)  
    
    ngdir_total = experi_common.ngdir_total;
    %%%% We obtain the multiple diffusion-encoding directions uniformly
    %%%% distributed in the unit sphere in 3 dimensions.
    [points_gdir,graddir_index,negii] = HARDI_PTS(ngdir_total);

   
    if (~isempty(experi_btpde))                      
        %%%% We solve the BTPDE in many diffusion directions.
        [SIG_BTPDE_cmpts_hardi,SIG_BTPDE_allcmpts_hardi,ctime_btpde_hardi] ...
            = SIG_BTPDE_HARDI(experi_btpde,mymesh,DIFF_cmpts,kappa_bdys,IC_cmpts,...
            points_gdir,graddir_index,negii);
        %%%% We plot the normalized signal in many diffusion directions.
        if (DO_PLOTS)    
            S0 = sum((IC_cmpts.*VOL_cmpts));         
            nb = size(experi_btpde.bvalues,2);                   
            for iexperi = 1:nexperi
                for ib = 1:nb                  
                    bv = experi_btpde.bvalues(iexperi,ib);
                    title_str = ['BTPDE. All Cmpts. ','Experi ',...
                        num2str(iexperi),', b= ',num2str(bv)];
                    PLOT_HARDI_PT(points_gdir,...
                        squeeze(real(SIG_BTPDE_allcmpts_hardi(:,iexperi,ib)))/S0,title_str);
                end
            end
        end
    end    
    %%%%% We solve the HADC model
    if (~isempty(experi_hadc))        
        %%%% We solve the HADC model in many diffusion directions.
        [ADC_HADC_cmpts_hardi,ADC_HADC_allcmpts_hardi,ctime_hadc_hardi] ...
            = HADC_HARDI(experi_hadc,mymesh,DIFF_cmpts,IC_cmpts,...
            points_gdir,graddir_index,negii);      
        %%%% We plot the normalized ADC (ADC/D0) in many diffusion directions.
        if(DO_PLOTS)
            ADC0 = sum((DIFF_cmpts.*VF_cmpts));
            for iexperi = 1:nexperi
                title_str = ['HADC. All Cmpts. ','Experi ',num2str(iexperi)];
                PLOT_HARDI_PT(points_gdir,...
                    squeeze(ADC_HADC_allcmpts_hardi(:,iexperi)/ADC0),title_str);
            end
        end
    end
end

