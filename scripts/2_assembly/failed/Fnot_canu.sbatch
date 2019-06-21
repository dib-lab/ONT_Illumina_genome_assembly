#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_canu
#SBATCH -J Fnot-canu
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_canu/canu-Fnot%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_canu/canu-Fnot%j.o
#SBATCH -t 336:00:00
#SBATCH -p LM
#SBATCH --mem=1000GB
#SBATCH --ntasks-per-node 14
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_canu

canu useGrid=false -p Fnot -d Fnot-canu genomeSize=1.1m gnuplotTested=true -nanopore-raw F_nottii.fastq

echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU,MaxDiskWrite,MaxDiskRead,NTasks,ConsumedEnergy' -P ${SLURM_JOB_ID}.batch | /home/ljcohen/parse_acc.py
