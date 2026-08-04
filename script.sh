#!/bin/bash

#SBATCH --job-name=BIODICAPT
#SBATCH --partition=std
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3 # number of chains to run (NOT number of K-folds)
#SBATCH --time=00-03:00:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

## Load software modules
module load R/4.6.1/gcc-openmpi
cd "$SLURM_SUBMIT_DIR"

## Stop if there is an error
set -e

## Clear outputs
rm -rf outputs/results
rm -rf outputs/figures

## Run scripts
Rscript -e 'source(".Rprofile"); source("scripts/3-STOC-OED.R")'

## Save datetime of run
cat > outputs/RUN.md <<EOF
Generated on: $(date '+%Y-%m-%d %H:%M:%S %Z')
EOF

## Finish gracefully
exit 0