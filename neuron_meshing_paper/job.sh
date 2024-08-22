#!/bin/bash

#SBATCH --job-name=neuron_meshing_paper
#SBATCH --output=neuron_meshing_paper/job_mf_%j.txt
#SBATCH --error=neuron_meshing_paper/error_mf_%j.txt

#SBATCH --time=12:00:00
#SBATCH --ntasks=40 ## number of nodes required
## load modules
module load matlab
## execution
# meshes=$(<neuron_meshing_paper/cells.txt);
# setup_file='setup_bvalues_experiments'
# tet="-pq1.2a1.0VCn";

# for mesh in $meshes; do
#      matlab -r "driver_segmented('$mesh','$setup_file','$tet','1.0')";       
#      echo "Completed $mesh";
# done;

# mesh=mesh_files/microglial/20-sn-1.CNG.ply
mesh=mesh_files/selected/1-2-1.CNG.ply

setup_file='setup_bvalues_experiments'
tet="-pq1.2a1.0VCn";

matlab -r "driver_cell('$mesh','$setup_file','$tet','3.0')";       
echo "Completed $mesh";


