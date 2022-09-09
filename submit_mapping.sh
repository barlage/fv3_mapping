#!/bin/bash

#BATCH --job-name=fv3_mapping
#SBATCH -t 00:10:00
#SBATCH -A bigmem
#SBATCH -A gsienkf
#SBATCH --qos=batch
#SBATCH -o mapping.out
#SBATCH -e mapping.out
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1

./create_fv3_mapping.exe



