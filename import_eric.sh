#!/bin/bash
#SBATCH --partition bigmem 
#SBATCH --spread-job
#SBATCH -n 1 -N1
#SBATCH --cpus-per-task 5
#SBATCH --mem 500G
#SBATCH --time 1-0:00:00
#SBATCH --job-name PapaBigmemJob
#SBATCH --output Papa%J.log
#SBATCH --mail-type=ALL
#SBATCH --mail-user=eric.song@yale.edu


echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------------'
echo -e '\n\n'


source activate R_merge_10x

mpirun R --slave -f eric2.r


closeCluster(cl)
mpi.quit()

echo '-------------------------------'
cd ${SLURM_SUBMIT_DIR}
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo '-------------------------------'
echo -e '\n\n'
q()

