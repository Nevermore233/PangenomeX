#!/bin/bash
#SBATCH -N 1                        # Request 1 node
#SBATCH -n 28                       # Request 28 tasks
#SBATCH --ntasks-per-node=28         # Specify 28 tasks per node
#SBATCH --output=%j.out              # Standard output will go to jobID.out
#SBATCH --partition=normal1         # Job partition
#SBATCH --error=%j.err               # Standard error will go to jobID.err

# ‘À–– Python Ω≈±æ
python3 tree2graph.py -nwk mynwk.nwk -npz data/bub_results.npz -cnv /data/home/std_13/000_pgcnv/data/zipcall-output/zipcaller_res_2025-04-01_16-25-47.cnv -o graph.gpickle






















