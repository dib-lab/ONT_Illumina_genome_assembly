#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_ONT
#SBATCH -J Fnot-nanofilt
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_ONT/nanofilt-Fnot%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_ONT/nanofilt-Fnot%j.o
#SBATCH -t 12:00:00
#SBATCH -p RM
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_ONT

gunzip -c F_nottii_trimmed.fq.gz | NanoFilt -q 5 | gzip > F_nottii_trimmed.Q5.fq.gz
