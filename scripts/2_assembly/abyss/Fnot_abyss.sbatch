#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_abyss
#SBATCH -J Fnot-abyss
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_abyss/abyss-Fnot%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_abyss/abyss-Fnot%j.o
#SBATCH -t 336:00:00
#SBATCH -p LM
#SBATCH --mem=1000GB
#SBATCH --ntasks-per-node 14
#SBATCH --cpus-per-task 1
source ~/.bashrc
source activate assembly
cd /pylon5/mc5phkp/ljcohen/kfish_abyss

SPECIES=Fnot
PROJECTDIR=$LOCAL/$SPECIES
mkdir $PROJECTDIR
cd $PROJECTDIR
cp /pylon5/bi5fpmp/ljcohen/kfish_Illumina/trimmed/Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L6_1.qc.fq.gz .
cp /pylon5/bi5fpmp/ljcohen/kfish_Illumina/trimmed/Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L6_2.qc.fq.gz .
cp /pylon5/bi5fpmp/ljcohen/kfish_Illumina/trimmed/Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L7_1.qc.fq.gz .
cp /pylon5/bi5fpmp/ljcohen/kfish_Illumina/trimmed/Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L7_2.qc.fq.gz .

abyss-pe k=51 name=Fcat_abyss in="Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L6_1.qc.fq.gz Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L6_2.qc.fq.gz Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L7_1.qc.fq.gz Fnotti_USPD16092508-N705-AK390_HV3JCCCXY_L7_2.qc.fq.gz" contigs

# Grab the assembly files and copy it to your storage
mkdir /pylon5/mc5phkp/ljcohen/kfish_abyss/Fnot
cp $PROJECTDIR/* /pylon5/mc5phkp/ljcohen/kfish_abyss/Fnot

# remove temp files on $LOCAL
rm -rf $PROJECTDIR

echo ${SLURM_JOB_NODELIST}       # Output Contents of the SLURM NODELIST

env | grep SLURM            # Print out values of the current jobs SLURM environment variables

scontrol show job ${SLURM_JOB_ID}     # Print out final statistics about resource uses before job exits

sstat --format 'JobID,MaxRSS,AveCPU,MaxDiskWrite,MaxDiskRead,NTasks,ConsumedEnergy' -P ${SLURM_JOB_ID}.batch | /home/ljcohen/parse_acc.py
