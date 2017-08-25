#!/bin/bash

#SBATCH --time=48:00:00

#SBATCH -n 1
#SBATCH --mem=8G

#SBATCH -J MIBI_5s_18
#SBATCH -o MIBI_5s_18.out
#SBATCH -e MIBI_5s_18.out

#SBATCH --dependency=afterany:10205876

matlab -r 'master_script; exit'

echo Finished execution at `date`
