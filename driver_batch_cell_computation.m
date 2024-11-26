%DRIVER_SPINDOCTOR Solve BTPDE, MF for many cells.
%   The user is advised to read the latest version
%   from \url{https://github.com/jingrebeccali/SpinDoctor}

clear
restoredefaultpath;

% Add SpinDoctor
addpath(genpath("src"));


%% Define inputs

% Get setup
setup_file='setup_batch';

meshes = ["mesh_files\selected\1-2-2.CNG.ply",
    "mesh_files\spindle\04b_spindle3aFI.ply",
    "mesh_files\Alzheimer_study\Ramified\ctrl_010319_13_826-2_1.ply"
];
% Different cells may need different refinement levels for the meshes.
% the default from the setup file can be overriden here/
tetgen_options= [
  "-pq1.2a0.5O9VCn",
  "-pq1.2a1.0O9VCn",
  "-pq1.2a0.5O9VCn"
];
ncells = length(meshes);

[~,cellnames,~] = fileparts(meshes);

%% Run simulations
% Run only once to obtain signals and finite element meshes
for i = 1:ncells
    mesh = meshes(i); tet_opt = tetgen_options(i);
    [results,femesh_cell,~,~]= run_simulations_neuron(mesh,setup_file,tet_opt);
end

%% Load simulations
addpath(genpath('drivers_postprocess'));
% Load signals and finite element meshes.
for i = 1:ncells
    mesh = meshes(i); tet_opt = tetgen_options(i);
    [results,femesh_cell,~,~]= load_simulations_neuron(mesh,setup_file,tet_opt);
    % Do any data analysis, visualisation or storage here.
    signals = results.mf_cell.signal/femesh_cell.total_volume;
    figure;
    subplot(1,2,1);
    plot_hardi_shells(results.setup,signals);
    title('Hardi plot')
    
    subplot(1,2,2);
    plot_femesh(femesh_cell,1,false);
    title('Finite element mesh');
    sgtitle(cellnames(i),'Interpreter','none');
end
