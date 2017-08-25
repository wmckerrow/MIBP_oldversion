#!/bin/bash

#SBATCH --time=48:00:00

#SBATCH -n 1
#SBATCH --mem=8G

#SBATCH -J conflicting_5s_18
#SBATCH -o conflicting_5s_18.out
#SBATCH -e conflicting_5s_18.out

#SBATCH --dependency=afterany:10205876

matlab -r "pg_eval; load([RNA_NAME '.mat']); prob_of_conflicting; exit"

echo Finished execution at `date`
