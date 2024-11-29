#!/bin/bash
#
#SBATCH --job-heatmap
#SBATCH --output=data.mat
#SBATCH --ntasks-per-node=40
#SBATCH --nodes=1
#SBATCH --time=30:00
#SBATCH -p extended-40core


module load module load matlab/2022a
matlab -nodisplay -nodesktop

cd /gpfs/home/armishra/cft
run heatmap.m
