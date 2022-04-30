#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH --partition=gpu
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each. 
#SBATCH -n 1
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=4:00:00
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err
source /sw/packages/python-tools/anaconda3/etc/profile.d/conda.sh
conda activate mamba
python init.py
signac view
python src/engines/gomc/project.py status
python src/engines/gomc/project.py submit
