#!/bin/bash
# Example with 28 cores for OpenMP
#
# Project/Account
#SBATCH --partition=gpu
# Request one node
#SBATCH -N 1
# Total number of cores, in this example it will 1 node with 1 core each. 
#SBATCH -n 5
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
# Runtime of this jobs is less then 12 hours.
#SBATCH --time=4:00:00
#SBATCH --mail-user=go2432@wayne.edu

#SBATCH -o output_%j.out

#SBATCH -e errors_%j.err
source /sw/packages/python-tools/anaconda3/etc/profile.d/conda.sh
export PYTHONPATH=$PYTHONPATH:/home/udpk5hr/reproducibility_study/reproducibility_project/:/home/udpk5hr/reproducibility_study
conda activate mamba
module load cuda
nvidia-smi
python src/engines/namd/project.py run

