#!/bin/bash

#SBATCH --job-name=TEST
#SBATCH --partition=std
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --time=00-06:00:00

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

## Load software modules
module load R/4.6.1/gcc-openmpi

## Stop if there is an error
set -e

## Clear outputs
rm -rf outputs/results
rm -rf outputs/figures

## Run scripts
Rscript scripts/3-STOC-OED.R

## Finish gracefully
exit 0