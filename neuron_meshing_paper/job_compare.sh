#!/bin/bash

#SBATCH --job-name=idefix_compare
#SBATCH --output=neuron_meshing_paper/output_compare_%j.txt
#SBATCH --error=neuron_meshing_paper/error_compare_%j.txt

#SBATCH --time=72:00:00
#SBATCH --ntasks=60 ## number of nodes required
## load modules
module load matlab
## execution
# meshes=$(<neuron_meshing_paper/cells.txt);
# meshes=mesh_files/selected/*.ply
setup_file='setup_comparison'
tet="-pq1.2a0.1VCn";
meshes=mesh_files/ultraliser/1-2-1-watertight.ply
for mesh in $meshes; do
     matlab -r "driver_btpde('$mesh','$setup_file','$tet','1')";       
     echo "Completed $mesh";
done;