#!/bin/bash

#SBATCH --job-name=idefix_cell_sim
#SBATCH --output=job_outputs/output_cell_sim_%j.txt
#SBATCH --error=job_errors/error_cell_sim_%j.txt

#SBATCH --time=2:00:00
#SBATCH --ntasks=10 ## number of nodes required
## load modules
module load matlab
## execution
meshes=$(<cells_human_separated.txt);
setup_file='setup_morez'
tet="-pq1.2a0.5VCn";

for mesh in $meshes; do
     matlab -r "driver_cell('$mesh','$setup_file','$tet','1.0')";       
     echo "Completed $mesh";
done;
