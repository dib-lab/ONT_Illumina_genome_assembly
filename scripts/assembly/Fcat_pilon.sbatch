#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat
#SBATCH -J Fcat.pilon
#SBATCH -t 168:00:00
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat/pilon-Fcat%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat/pilon-Fcat%j.o
#SBATCH -p LM
#SBATCH --mem=2000GB
source ~/.bashrc
source activate ONT_pilon 
module load samtools/1.9
samtools index Fcat_pilon.sorted.bam
java -jar -Xmx1500G /pylon5/bi5fpmp/ljcohen/miniconda3/envs/ONT_pilon/share/pilon-1.23-0/pilon-1.23.jar --genome Fcat.ctg.fa --bam Fcat_pilon.sorted.bam --threads 32

echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU,MaxDiskWrite,MaxDiskRead,NTasks,ConsumedEnergy' -P ${SLURM_JOB_ID}.batch | /home/ljcohen/parse_acc.py
