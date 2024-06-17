addpath(genpath("src"));


%% Define inputs

% Get setup
addpath(genpath("setups"));
setup_dendrite_soma;

%% Prepare simulation
[setup, femesh, ~, ~]  = prepare_simulation(setup);

%% Perform small experiments
% Short time approximation (STA) of the ADC
% Free diffusion signal
free = compute_free_diffusion(setup.gradient.bvalues, setup.pde.diffusivity, ...
    femesh.volumes, setup.pde.initial_density);
[~,cellname,~] =fileparts(setup.name); 
savepath_cell = sprintf("saved_simul/%s/cell_D_%f",cellname,setup.pde.diffusivity_out(1,1));

tic
% Perform BTPDE experiments
if isfield(setup, "btpde")
    % Solve BTPDE
    btpde = solve_btpde(femesh, setup,savepath_cell,false);
end
toc

%%
fprintf("Computing soma signal for %s",cellname);
file = sprintf("mesh_files/alphaSwc_meshes_small/%s_ply_dir/%s_no_ecs_tet%s_mesh.1",cellname,cellname,setup.geometry.tetgen_options);
neighbours =read_tetgen_neigh(file);
swc_file = fullfile("swc_files",cellname+".swc");
swc = Swc(swc_file);
femesh_soma = separate_soma(femesh,swc,neighbours);
savepath_soma  = sprintf("saved_simul/%s/soma_D%f",cellname,setup.pde.diffusivity_out(1,1));
setup.btpde.rerun = false;
tic
% Perform BTPDE experiments
if isfield(setup, "btpde")
    % Solve BTPDE
    btpde_soma = solve_btpde(femesh_soma, setup,savepath_soma,false);
end
toc



%%
disp("Separating out dendrite branches")
dendrites = find_dendrites_femesh(femesh,femesh_soma,neighbours);

disp("Computing dendrite signals");
for iden = 1:length(dendrites)
    femesh_dendrite = dendrites{iden};
    savepath_dendrite = sprintf("saved_simul/%s/dendrite_%d_D%f",cellname,iden,setup.pde.diffusivity_out(1,1));   
    tic
    % Perform BTPDE experiments
    if isfield(setup, "btpde")
        % Solve BTPDE
        btpde = solve_btpde(femesh_dendrite, setup,savepath_dendrite,false);
    end
    toc
end




