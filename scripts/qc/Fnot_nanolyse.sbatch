#!/bin/bash -l
#SBATCH -D /pylon5/mc5phkp/ljcohen/kfish_ONT/
#SBATCH -J Fnot-nanoly
#SBATCH -o /pylon5/mc5phkp/ljcohen/kfish_ONT/nanoly-Fcat%j.o
#SBATCH -e /pylon5/mc5phkp/ljcohen/kfish_ONT/nanoly-Fcat%j.o
#SBATCH -t 12:00:00
#SBATCH -p RM
#SBATCH --ntasks-per-node 4
#SBATCH --cpus-per-task 1

source ~/.bashrc
source activate ONT

cd /pylon5/mc5phkp/ljcohen/kfish_ONT

gunzip -c F_nottii.fq.gz | NanoLyse --reference ../ref/DNA_CS.fasta | gzip > F_nottii_woLambda.fq.gz
