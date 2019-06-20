#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_wtdbg
#SBATCH -J Foli-wtdbg
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_wtdbg/wtdbg-Foli%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_wtdbg/wtdbg-Foli%j.o
#SBATCH -t 48:00:00
#SBATCH -p LM
#SBATCH --mem=200GB
#SBATCH --ntasks-per-node 14
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_wtdbg

wtdbg2 -x ont -g1.1g -t 12 -i F_olivaceus_trimmed.Q5.fq.gz -fo Foli
wtpoa-cns -t 12 -i Foli.ctg.lay.gz -fo Foli.ctg.fa


echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU,MaxDiskWrite,MaxDiskRead,NTasks,ConsumedEnergy' -P ${SLURM_JOB_ID}.batch | /home/ljcohen/parse_acc.py
