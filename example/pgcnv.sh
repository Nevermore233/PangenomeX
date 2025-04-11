#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --partition=p100
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# ‘À–– Python Ω≈±æ
python3 pgcnv.py -config parameter_cfg_20250225.config -k  /data/home/std_13/000_pgcnv/graph.gpickle -o data/pgcnv_output






















