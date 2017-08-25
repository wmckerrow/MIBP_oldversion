#!/bin/bash

#SBATCH --time=48:00:00

#SBATCH -n 1
#SBATCH --mem=8G

#SBATCH -J disp_5s_18
#SBATCH -o disp_5s_18.out
#SBATCH -e disp_5s_18.out

#SBATCH --dependency=afterany:10205876

matlab -r 'redisplay; exit'

echo Finished execution at `date`
