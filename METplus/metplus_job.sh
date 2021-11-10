#!/bin/bash
#SBATCH --account=ovp
#SBATCH --job=metplus_job
#SBATCH --output=/scratch2/NCEPDEV/ovp/Marcel.Caron/METplus/output/jobs_clo/metplus_job.%j.out
#SBATCH --nodes=10
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --partition=hera

set -x

conf_base="/scratch2/NCEPDEV/stmp1/Marcel.Caron/scripts/s2s_scripts/METplus/"
metplus_conf_fname="CFS_PCP_grid_stat_metplus.conf"
hera_conf=${conf_base}/hera.conf
srun --ntasks 1 --nodes 1 --exclusive master_metplus.py \
   -c ${conf_base}/${metplus_conf_fname} \
   -c ${hera_conf} \

wait
sleep 5
