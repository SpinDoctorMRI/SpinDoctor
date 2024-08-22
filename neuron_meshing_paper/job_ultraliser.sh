#!/bin/bash

#SBATCH --job-name=neuron_meshing_paper
#SBATCH --output=neuron_meshing_paper/job_ultraliser_%j.txt
#SBATCH --error=neuron_meshing_paper/error_ultraliser_%j.txt

#SBATCH --time=24:00:00
#SBATCH --ntasks=40 ## number of nodes required
## load modules
module load matlab
## execution
mesh=mesh_files/ultraliser/1-2-1-watertight.ply 
setup_file='setup_comparison'
tet="-pq1.2a0.1VCn";

matlab -r "driver_cell('$mesh','$setup_file','$tet','3.0')";       
echo "Completed $mesh";

# mesh=mesh_files/ultraliser/20-sn-1-watertight.ply
# matlab -r "driver_cell('$mesh','$setup_file','$tet','1.0')";       
# echo "Completed $mesh";


# mesh=mesh_files/selected/02a_pyramidal2aFI.CNG.ply
# setup_file='setup_bvalues_experiments'
# tet="-pq1.2a1.0VCn";

# matlab -r "driver_segmented('$mesh','$setup_file','$tet','1.0')";       
# echo "Completed $mesh";


