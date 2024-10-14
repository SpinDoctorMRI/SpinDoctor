cd ..
addpath(genpath('src'))
addpath(genpath('setups'))
setup_file='setup_pgse_ref_sol';
tet="-pq1.2a0.5O9VCn";
tet_btpde="-pq1.2a0.1O9VCn";

fid = fopen("cells_human.txt",'r');
tline = fgetl(fid);
i=1;
clear meshes;
while ischar(tline)
    meshes(i) = string(tline);
    tline = fgetl(fid);
    i = i + 1;
end
ncells = length(meshes)
clear types;

for i = 1:ncells
    file = meshes(i);
    sep_file = split(file,"/");
    types(i) = sep_file(3);
end

clear errors; icell = 0;
clear volumes; clear volumes_ult;
clear meshes_tested;
for i =1:ncells
    mesh = meshes(i); type = types(i)
    [~,cell,~] = fileparts(mesh); 
    mesh_ult = sprintf("mesh_files/ultraliser_modified/%s/%s_um.ply",type,cell);
    if isfile(mesh_ult)
    results_neuron= load_mf_cell(mesh,setup_file,tet,'1.0');
    results_neuron_ultraliser = load_mf_cell(mesh_ult,setup_file,tet,'1.0');
    bvals = results_neuron.setup.gradient.bvalues;
    save(sprintf('neuron_meshing_paper/output_data/%s.mat',cell),'results_neuron','results_neuron_ultraliser','bvals');
    signal = real(results_neuron.mf.signal/results_neuron.total_volume);
    signal_ultraliser = real(results_neuron_ultraliser.mf.signal/results_neuron_ultraliser.total_volume);
    n = length(bvals);
    error = abs((signal- signal_ultraliser)./signal);
    icell = icell + 1; errors(icell,:) = error; meshes_tested(icell) = cell;
    disp(error);
    volumes(icell) = results_neuron.total_volume; volumes_ult(icell) = results_neuron_ultraliser.total_volume;
    end
end
rel_vol_errors = abs(volumes-volumes_ult)./volumes
