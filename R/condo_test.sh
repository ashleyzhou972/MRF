#!/bin/bash

# Copy/paste this job script into a text file and submit with the command:
#    sbatch thefilename
# job standard output will go to the file slurm-%j.out (where %j is the job ID)

#SBATCH --time=50:00:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1   # number of nodes
#SBATCH --ntasks-per-node=16   # 16 processor core(s) per node 
#SBATCH --mem=32G   # maximum memory per node
#SBATCH --job-name="tuning_test"
#SBATCH --mail-user=nzhou@iastate.edu   # email address
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
time /work/friedberg_lab/nzhou/configure/R_make/bin/Rscript tuning_real_3_condo.R  20180731 20
# First try two iterations to 
#	1 calculate eta domain
#	2 estimate time per iteration
# Then try two iteration to 
#	1 tune variance for jump proportion
#	2 further estimate time per iteration
