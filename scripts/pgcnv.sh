#!/bin/sh
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --partition=p100
#SBATCH --output=%j.out
#SBATCH --error=%j.err

# ‘À–– Python Ω≈±æ
python3 pgcnv.py -config my.config -k graph.gpickle -o data/pgcnv_output






















