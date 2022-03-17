#!/bin/bash
#SBATCH --account=ovp
#SBATCH --job=hght_anom
#SBATCH --output=/scratch2/NCEPDEV/stmp1/Shannon.Shields/scripts/s2s/cases/sbatch_out/python_job.%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --partition=hera

set -x

# Summary: Run any python script as a batch job.
# If it doesn't work on the first try, sometimes trying again will work
# Still may not work for every python script.
# To use:
# 1) Change line 4 to desired directory for log output
# 2) Write the location of python script that you want to run (line 20)
# 3) Write the name of the desired python script in the srun command (line 21)
# 4) run this file with sbatch on the command line
script_dir="/scratch2/NCEPDEV/stmp1/Shannon.Shields/scripts/s2s/repo_scripts"
srun python ${script_dir}/height_anoms.py
sleep 5
