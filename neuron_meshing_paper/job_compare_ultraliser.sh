#!/bin/bash

#SBATCH --job-name=idefix_compare
#SBATCH --output=neuron_meshing_paper/output_compare_%j.txt
#SBATCH --error=neuron_meshing_paper/error_compare_%j.txt

#SBATCH --time=12:00:00
#SBATCH --ntasks=20 ## number of nodes required
## load modules
module load matlab;
## execution

cd neuron_meshing_paper;
matlab -batch "compare_ultraliser_signals";