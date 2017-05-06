#!/bin/bash
#SBATCH --job-name Rjob      # Set a name for your job. This is especially useful if you have multiple jobs qeued
#SBATCH --partition student     # Slurm partition to use
#SBATCH --ntasks 1          # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 0-12:00        # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=1024     # Memory limit for each tasks (in MB)
#SBATCH -o Rjob.out    # File to which STDOUT will be written
#SBATCH -e Rjob.err    # File to which STDERR will be written
#SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=yxcheng@uw.edu # Email to which notifications will be sent

R CMD BATCH bayeslogisticreg-startcode-HW4.R
