#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=52
#SBATCH --mem=100000
#SBATCH --time=12:00:00

#SBATCH --mail-user=neil.g.maclaren@gmail.com
#SBATCH --mail-type=ALL
## #SBATCH --job-name=fMRI
## #SBATCH --output=./output/fMRI.out
#SBATCH --partition=general-compute
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

# for timing purposes, constrain the hardware
#SBATCH --constraint=CPU-Gold-6448Y

module load gcc openmpi r

dynamics=$1

Rscript ../sims/collect-fmri.R --dynamics=$dynamics
