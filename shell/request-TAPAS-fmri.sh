#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4000
#SBATCH --time=06:00:00

#SBATCH --mail-user=neil.g.maclaren@gmail.com
#SBATCH --mail-type=NONE
#SBATCH --qos=general-compute
#SBATCH --cluster=ub-hpc

#SBATCH --constraint=CPU-Gold-6448Y

module load matlab/2023b

# infile=\'/projects/academic/naokimas/neil/brains-ns50/101107.txt\'

infile=$1
echo ${infile}

matlab -batch "TAPASfmri($infile)"
