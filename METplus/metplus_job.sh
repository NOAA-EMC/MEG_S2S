#!/bin/bash
#SBATCH --account=ovp
#SBATCH --job=metplus_job
#SBATCH --output=/scratch2/NCEPDEV/stmp1/Shannon.Shields/scripts/s2s/cases/met_out/metplus_job.%j.out
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --partition=hera

set -x

conf_base="/scratch2/NCEPDEV/stmp1/Shannon.Shields/scripts/s2s/repo_scripts/METplus/"
metplus_conf_fname="CFS_TMP_grid_stat_metplus.conf"
hera_conf=${conf_base}/hera.conf
srun --ntasks 1 --nodes 1 --exclusive run_metplus.py \
   -c ${conf_base}/${metplus_conf_fname} \
   -c ${hera_conf} \

wait
sleep 5
