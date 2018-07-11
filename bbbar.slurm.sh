#!/bin/bash
#SBATCH --job-name=jdb12_pythia_bbbar
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-01:20:00
#SBATCH --partition=scavenge
#SBATCH --export=ALL

cd $SHARED_SCRATCH/jdb12/pythia_bbbar/
srun /home/jdb12/work/PythiaBBbar/run.sh
# srun /home/jdb12/work/Pythia6TreeMaker/GENERATOR 4 1000000 $SLURM_ARRAY_TASK_ID
