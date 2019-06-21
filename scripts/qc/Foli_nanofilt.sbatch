#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_ONT
#SBATCH -J Foli-nanofilt
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_ONT/nanofilt-Foli%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_ONT/nanofilt-Foli%j.o
#SBATCH -t 12:00:00
#SBATCH -p RM
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_ONT

gunzip -c F_olivaceus_trimmed.fq.gz | NanoFilt -q 5 | gzip > F_olivaceus_trimmed.Q5.fq.gz
