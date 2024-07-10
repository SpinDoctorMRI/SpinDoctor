#!/bin/bash

#SBATCH --job-name=idefix_cell_sim
#SBATCH --output=job_outputs/output_cell_sim_%j.txt
#SBATCH --error=job_errors/error_cell_sim_%j.txt

#SBATCH --time=96:00:00
#SBATCH --ntasks=64 ## number of nodes required
## load modules
module load matlab
## execution
meshes=$(<cells_human.txt);
tetgen_options="-pq1.2a0.1VCn";
setup_file='setup_morez_ref_sol';

for mesh in $meshes; do
    matlab -r "driver_btpde('$mesh','$setup_file','$tetgen_options','1')";
    echo "Completed $mesh";
done;

