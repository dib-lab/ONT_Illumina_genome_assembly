#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat
#SBATCH -J Fcat_bwa
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat/bwa-Fcat%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_pilon/Fcat/bwa-Fcat%j.o
#SBATCH -t 36:00:00
#SBATCH -p LM
#SBATCH --mem=500GB
module load samtools/1.9 bwa/0.7.13 
bwa index Fcat.ctg.fa
bwa mem -t 32 Fcat.ctg.fa Fcat_L6L7_1.qc.fq.gz Fcat_L6L7_2.qc.fq.gz | samtools view -Sb - | samtools sort - -o Fcat_pilon.sorted.bam


echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU,MaxDiskWrite,MaxDiskRead,NTasks,ConsumedEnergy' -P ${SLURM_JOB_ID}.batch | /home/ljcohen/parse_acc.py
