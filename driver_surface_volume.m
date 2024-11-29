%DRIVER_SURFACE_VOLUME Tries to estimate surface volume ratio from signals.
%
%   Provides a template to demonstrate how camino scheme sequences can be
%   integrated into SpinDoctor
%
%   The user is advised to read the latest version
%   from \url{https://github.com/jingrebeccali/SpinDoctor}

% Add SpinDoctor
addpath(genpath("src"));addpath(genpath('setups'));
addpath(genpath('src'));

% Read mesh paths
meshes = readlines("cells_human.txt");
ncells = length(meshes);

%% Set experiment parameters

tetgen_options = "-pq1.2a0.05O9VCn";
setup_file = "setup_surface_volume_ratio";


%% Run simulations fpr PGSE sequences
% Run only once to obtain signals and finite element meshes
for i = 1:ncells
    mesh = meshes(i); 
    [results,femesh_cell,~,~]= run_simulations_neuron(mesh,setup_file,tetgen_options);
end

%% Load simulations for PGSE results, Classical mitra formula

nsequence = 3;
for iseq = 1:nsequence
figure(iseq);
hold on; grid on;
xlabel("$S/V \, (\mu m^{-1})$",'Interpreter','latex','FontSize',20)
ylabel("$\frac{-9(D_{msd}/D_0  - 1)\sqrt{\pi}}{4\sqrt{D_0t}} \, (\mu m^{-1})$",'Interpreter','latex','FontSize',20)
title(sprintf("Sequence %d",iseq));
end

for i = 1:ncells
    mesh = meshes(i); 
    [results,femesh,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
    V = femesh.total_volume; S= femesh.total_area;
    bvalues = results.setup.gradient.bvalues(:,1);
    sequences =  results.setup.gradient.sequences;

    for iseq = 1:length(sequences)
        fig = figure(iseq); hold on;
        seq = sequences{iseq};
        D0 = results.setup.pde.mean_diffusivity;
        signal = real(results.mf_cell.signal(:,:,iseq,:))/V;
        signal_allcmpts = real(results.mf_cell.signal_allcmpts(:,iseq,:))/V;
    
        % MITRA formula
        adc_fitting = fit_signal(signal, signal_allcmpts, bvalues);
        t = seq.Delta;
        D_msd = adc_fitting.adc;
        SV_ratio_dmri = -(D_msd/D0  - 1)*(3*3*sqrt(pi)/(4*sqrt(D0*t)));
        
        fprintf("Sequence: %s\nReal surface-volume ratio = %f, Predicted surface-volume ratio = %f\n",seq,S/V,SV_ratio_dmri);
    
        % Plotting
        is_amoeboid = contains(mesh,'Amoeboid','IgnoreCase',true);
        is_ramified = contains(mesh,'Ramified','IgnoreCase',true);
    
        if is_amoeboid
        
            if label_amoeboid
                scatter(S/V,SV_ratio_dmri,'r',"filled","DisplayName","Amoeboid");
                label_amoeboid= false;
            else
                scatter(S/V,SV_ratio_dmri,'r',"filled");
            end
        elseif is_ramified
        
            if label_ramified
                    scatter(S/V,SV_ratio_dmri,'b',"filled","DisplayName","Ramified");
                    label_ramified= false;
            else
                scatter(S/V,SV_ratio_dmri,'b',"filled");
            end
        end
    end
end
axis equal;
legend;

%% Load simulations for PGSE results, Updated mitra formula

nsequence = 3;
for iseq = (1:nsequence) + nsequence
figure(iseq);
hold on; grid on;
xlabel("$S/V \, (\mu m^{-1})$",'Interpreter','latex','FontSize',20)
ylabel("$\frac{-9(D_{msd}/D_0  - 1)\sqrt{\pi}}{4\sqrt{D_0t}\eta} \, (\mu m^{-1})$",'Interpreter','latex','FontSize',20)
title(sprintf("Sequence %d",iseq));
end

for i = 1:ncells
    mesh = meshes(i); 
    [results,femesh,~,~]= load_simulations_neuron(mesh,setup_file,tetgen_options);
    V = femesh.total_volume; S= femesh.total_area;
    bvalues = results.setup.gradient.bvalues(:,1);
    sequences =  results.setup.gradient.sequences;
    
    S3 = get_mitra_S3(femesh);

    for iseq = 1:length(sequences)
        fig = figure(iseq+nsequence); hold on;
        seq = sequences{iseq};
        D0 = results.setup.pde.mean_diffusivity;
        signal = real(results.mf_cell.signal(:,:,iseq,:))/V;
        signal_allcmpts = real(results.mf_cell.signal_allcmpts(:,iseq,:))/V;
    
        % UPDATED MITRA formula
        adc_fitting = fit_signal(signal, signal_allcmpts, bvalues);
        t = seq.Delta;
        D_msd = adc_fitting.adc;

        T3 = seq.get_T_m(3,results.setup.gamma);
        eta = trace(S3*T3);
        SV_ratio_dmri = -(D_msd/D0  - 1)*(3*sqrt(pi)/(4*sqrt(D0*T)*eta));
        
        fprintf("Sequence: %s\nReal surface-volume ratio = %f, Predicted surface-volume ratio = %f\n",seq,S/V,SV_ratio_dmri);
    
        % Plotting
        is_amoeboid = contains(mesh,'Amoeboid','IgnoreCase',true);
        is_ramified = contains(mesh,'Ramified','IgnoreCase',true);
    
        if is_amoeboid
        
            if label_amoeboid
                scatter(S/V,SV_ratio_dmri,'r',"filled","DisplayName","Amoeboid");
                label_amoeboid= false;
            else
                scatter(S/V,SV_ratio_dmri,'r',"filled");
            end
        elseif is_ramified
        
            if label_ramified
                    scatter(S/V,SV_ratio_dmri,'b',"filled","DisplayName","Ramified");
                    label_ramified= false;
            else
                scatter(S/V,SV_ratio_dmri,'b',"filled");
            end
        end
    end
end
axis equal;
legend;
